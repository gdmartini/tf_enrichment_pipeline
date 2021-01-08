% DESCRIPTION
% Script to conduct HMM inference
%
% ARGUMENTS
% project: master ID variable 
%
% wInf: memory used for inference
%
% KInf: number of states used for inference
%
% OPTIONS
% dropboxFolder: Path to data folder where you wish to save
%                pipeline-generated data sets and figures. If this
%                var is not specified, output will be saved one level
%                above git repo in the folder structure
%
% controlProject: specifies a project to use as an external control
%
% OUTPUT: hmm_input_output, structure containing vectors of protein and MS2
% intensities, along with corresponding HMM-decoded activity trajectories

function main06_incorporate_hmm_results(projectName,varargin)

close all
addpath(genpath('utilities'))

makeLongFormSet = 0;
myCluster = parcluster('local');
nWorkersMax = ceil(myCluster.NumWorkers/2);
bootstrap_flag = 0;
n_bootstraps = 1;

% check for optional inputs
for i = 1:(numel(varargin)-1)  
    if i ~= numel(varargin)        
        eval([varargin{i} '=varargin{i+1};']);                
    end    
end

% get path to results
if ~exist('resultsRoot','var')
    liveProject = LiveEnrichmentProject(projectName);
    resultsRoot = [liveProject.dataPath filesep];
    resultsDir = [resultsRoot 'cpHMM_results' filesep];
else
    resultsRoot = [resultsRoot filesep projectName filesep];
    resultsDir = [resultsRoot filesep 'cpHMM_results' filesep];
end

% get list of all inference subdirectories. By default, we'll generate
% summaries for all non-empty inference sub-directory
infDirList = dir([resultsDir 'w*']); % get list of all inference subdirectories. By default, we'll generate
   
for trace_i = 1:length(infDirList)
  
    % get path to results
    resultsPath = [infDirList(trace_i).folder filesep];
    
    % load inference options
    load([resultsPath infDirList(trace_i).name filesep 'inferenceOptions.mat'])
    if length(fieldnames(inferenceOptions))==1
        inferenceOptions = inferenceOptions.inferenceOptions;
    end
    
    % extract key inference parameters     
    nStates = inferenceOptions.nStates;
    nSteps = inferenceOptions.nSteps;
    alpha = inferenceOptions.alpha;
    Tres = inferenceOptions.Tres;
    timeGrid = 0:Tres:60*60;
    maxDT = 1.2*Tres; % maximum distance from observed data point
    
    % load compiled inference files
    load([resultsPath 'compiledResults_' infDirList(trace_i).name '.mat'],'compiledResults');
    
    % load corresponding trace structure
    if ~inferenceOptions.ProteinBinFlag
        load([resultsRoot filesep 'spot_struct.mat'])
        analysis_traces = spot_struct;
        clear spot_struct;
    else
        load([resultsRoot filesep 'spot_struct_protein.mat'])
        analysis_traces = spot_struct_protein;
        clear spot_struct_protein;
    end
               

    % check for existence of fit structure
    trace_fit_flag = 1;
    inf_files = dir([resultsPath infDirList(trace_i).name filesep 'hmm_results*.mat']);
    if exist([resultsPath 'singleTraceFits_' infDirList(trace_i).name '.mat']) > 0
        fit_props = dir([resultsPath 'singleTraceFits_' infDirList(trace_i).name '.mat']);
        fit_date = datenum(fit_props(1).date);
        hmm_date = datenum(inf_files(1).date);
        if fit_date > hmm_date
            trace_fit_flag = 0;
        end
    end
    
    % get list of particle IDs
    trace_particle_index = [analysis_traces.particleID];
    
    % perform viterbi trace decoding if necessary
    if trace_fit_flag    
        singleTraceFits = performSingleTraceFits(compiledResults, inferenceOptions, bootstrap_flag, ...
                  analysis_traces, trace_particle_index, nWorkersMax, resultsPath, infDirList(trace_i).name);        
    else
        load([resultsPath 'singleTraceFits_' infDirList(trace_i).name '.mat'],'singleTraceFits')
    end
  
    % generate longform dataset
    if false%makeLongFormSet        
        resultsTable = generateLongFormTable(analysis_traces, timeGrid, singleTraceFits, resultsPath, infDirList(trace_i).name);
    end

end

% generate alpha kernel for estimating predicted HMM fluroescence
alpha_kernel = NaN(1,nSteps);
for i = 1:nSteps
    if i < alpha
        alpha_kernel(i) = ((i-1) / alpha  + .5 * 1/alpha )*Tres;
    elseif i > alpha && (i-1) < alpha
        alpha_kernel(i) = Tres*(1 - .5*(alpha-i+1)*(1-(i-1)/alpha));
    else
        alpha_kernel(i) = Tres;
    end
end

%%% now combine with protein data and raw traces
disp('building input/output dataset...')

% first combine HMM and raw fluoresence trace data points
hmm_input_output = addTraceAnalysisFields(analysis_traces,singleTraceFits);

% next, incorporate local enrichment fields (if they exist)
if exist([resultsRoot filesep 'spot_struct_protein.mat'],'file')
    load([resultsRoot filesep 'spot_struct_protein.mat'],'spot_struct_protein')
    fitParticleVec = [singleTraceFits.particleID];
    traceParticleVec = [spot_struct_protein.particleID];
    
    for trace_i = 1:length(hmm_input_output)
          
        % basic indexing info
        ParticleID = hmm_input_output(trace_i).particleID(1);   
        fitIndex = find(fitParticleVec==ParticleID);
        traceIndex = find(traceParticleVec==ParticleID);
  
        % extract time vectors
        timeFrom = spot_struct_protein(traceIndex).time;
        timeTo = hmm_input_output(trace_i).time;
        
        % interpolate protein vectors to match cpHMM spacing
        interpFields = {'spot_protein_vec','serial_null_protein_vec', 'edge_null_protein_vec', 'nuclear_protein_vec',...
          'spot_mcp_vec'};
        
        for field_i = 1:length(interpFields)
            vecRaw = spot_struct_protein(traceIndex).(interpFields{field_i});
            nanFT = ~isnan(vecRaw);
            if sum(nanFT) > 1
                hmm_input_output(trace_i).(interpFields{field_i}) = interp1(timeFrom(nanFT), vecRaw(nanFT), timeTo);
            elseif sum(nanFT) == 1
                vecTo = NaN(size(timeTo));
                [~,mi] = min(abs(timeTo-timeFrom(nanFT)));
                vecTo(mi) = vecRaw(nanFT);
                hmm_input_output(trace_i).(interpFields{field_i}) = vecTo;
            else
                hmm_input_output(trace_i).(interpFields{field_i}) = NaN(size(timeTo));
            end
        end
                                             
        % generate flag var indicating interpolated obs that are too far from 
        % true points
        input_times = timeFrom(~isnan(spot_struct_protein(trace_i).spot_protein_vec));
        dt_vec_gap = NaN(size(hmm_input_output(trace_i).time));
        if ~isempty(input_times)
            for t = 1:numel(dt_vec_gap)
                dt_vec_gap(t) = min(abs(hmm_input_output(trace_i).time(t)-input_times));   
            end
        end
        hmm_input_output(trace_i).dt_filter_gap = dt_vec_gap > maxDT;            
        % record general info for later use        
        hmm_input_output(trace_i).Tres = Tres;
        hmm_input_output(trace_i).maxDT = maxDT;
        hmm_input_output(trace_i).InferenceID = fitIndex;            
        hmm_input_output(trace_i).rawSpotStructIndex = traceIndex;            
        
    end
end  

% save results
save([resultsDir 'hmm_input_output.mat'],'hmm_input_output')
