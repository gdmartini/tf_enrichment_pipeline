function cpHMMInferenceGrouped(projectName,InputDataPath,OutputDataPath)
% Script to call primary cpHMM wrapper function

clear
close all
warning('off','all') %Shut off Warnings
modelPath = './utilities';
addpath(modelPath);

% get inference options
inferenceOptions = determineInferenceOptions;

% project identifier
% projectName = '2xDl-Ven_hbP2P-mCh';
% default path to model scripts

% InputDataPath = ['../../dat/tf_enrichment/' projectName '/'];
 
%%%%%%%%%%%%%%



if inferenceOptions.ProteinBinFlag && inferenceOptions.savioFlag
    nBoots = 1; % will run multiple instances on savio
else  
    nBoots = 5;
end

% DataPath = ['S:\Nick\Dropbox\ProcessedEnrichmentData\' project '\'];%


% check that we have proper fields
if ~inferenceOptions.dpBootstrap
    warning('Bootstrap option not selected. Setting nBoots to 1')
    nBoots = 1;
end

%% %%%%%%%%%% Set fraction of gene length comprised of MS2 cassette %%%%%%%

if contains(projectName,'hbP2P')
    alphaFrac = 1275 / 4670;
elseif contains(projectName,'snaBAC')
    alphaFrac = 1302 / 6444;
end
alpha = alphaFrac*inferenceOptions.nSteps;

%% %%%%%%%%%%%%%%%%%%%%%% Load trace data set %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~inferenceOptions.ProteinBinFlag
  load([InputDataPath '/nucleus_struct.mat'],'nucleus_struct') % load data
  analysis_traces = nucleus_struct;
  clear nucleus_struct
else
  load([InputDataPath '/spot_struct_protein.mat'],'spot_struct_protein') % load data
  analysis_traces = spot_struct_protein;
  clear spot_struct_protein
end

%% %%%%%%%%%%%%%%%%%%%%%% Generate output directory %%%%%%%%%%%%%%%%%%%%%%%
% Set write path (inference results are now written to external directory)
if inferenceOptions.fluo3DFlag
    fluoSuffix = 'f3D';
else
    fluoSuffix = 'f2D';
end

if inferenceOptions.ProteinBinFlag
  outSuffix =  ['hmm_inference_protein' filesep 'w' num2str(inferenceOptions.nSteps) '_inferenceOptions.nStates' num2str(inferenceOptions.nStates) '_' fluoSuffix filesep]; 
else
  outSuffix =  ['hmm_inference_mf' filesep 'w' num2str(inferenceOptions.nSteps) '_inferenceOptions.nStates' num2str(inferenceOptions.nStates) '_' fluoSuffix filesep];   
end

% set write path
outPrefix = [OutputDataPath filesep projectName filesep]; %hmmm_data/inference_out/';
outDir = [outPrefix outSuffix];
mkdir(outDir);

%% %%%%%%%%%%%%%%%%%%%%%% Process/filter trace data %%%%%%%%%%%%%%%%%%%%%%%
trace_struct_filtered = filterTraces(inferenceOptions,analysis_traces);

% define iteration wrapper
iter_index = 1;
iter_ref_list = ones(size(trace_struct_filtered));
if inferenceOptions.ProteinBinFlag    
    iter_ref_list = [trace_struct_filtered.mf_protein_bin];
    iter_index = unique(iter_ref_list);
end

%% %%%%%%%%%%%%%%%%%%%%%% Conduct cpHMM Inference %%%%%%%%%%%%%%%%%%%%%%%%%
rng('shuffle'); % ensure we don't repeat bootstrap samples across different replicates

% iterate through designated groups
for t = 1:length(iter_index)
    iter_filter = iter_ref_list == iter_index(t);
    
    for b = 1:nBoots
      
        iter_start = now;
        local_struct = struct;    
        output = struct;                
        
        % Extract subset of traces relevant to this subgroup       
        inference_set = trace_struct_filtered(iter_filter);                
        set_size = length([inference_set.fluo]);  
        
        skip_flag = 0;
        if set_size < inferenceOptions.minDPperInf                    
            skip_flag = 1;                    
            warning('Too few data points. Skipping')                                    
        else 
            sample_index = 1:length(inference_set);
            
            %% take bootstrap sample
            ndp = 0;    
            sample_ids = [];                    
            
            %Reset bootstrap size to be on order of set size for small bins            
            inferenceOptions.SampleSize = min([inferenceOptions.SampleSize ceil(set_size/100)*100]);
            
            % randomly draw traces
            while ndp < inferenceOptions.SampleSize
                tr_id = randsample(sample_index,1);
                sample_ids = [sample_ids tr_id];
                ndp = ndp + length(inference_set(tr_id).time);
            end
            
            % add them to data cells
            fluo_data = cell([length(sample_ids), 1]);    
            time_data = cell([length(sample_ids), 1]);    
            sample_particles = [inference_set(sample_ids).ParticleID];
            for tr = 1:length(sample_ids)
                fluo_data{tr} = inference_set(sample_ids(tr)).fluo;                    
                time_data{tr} = inference_set(sample_ids(tr)).time;                    
            end            
           
            %% Random initialization of model parameters
            param_init = initialize_random (inferenceOptions.nStates, inferenceOptions.nSteps, fluo_data);
            % Approximate inference assuming iid data for param initialization                
            local_iid_out = local_em_iid_reduced_memory(fluo_data, param_init.v, ...
                                param_init.noise, inferenceOptions.nStates, inferenceOptions.nSteps, alpha, 500, 1e-4);
                              
            noise_iid = 1/sqrt(exp(local_iid_out.lambda_log));
            v_iid = exp(local_iid_out.v_logs);  
            
            %% create parallel pool if one does not already exist
            p = gcp('nocreate');
            if isempty(p)
                parpool(inferenceOptions.maxWorkers); %6 is the number of cores the Garcia lab server can reasonably handle per user.
            elseif p.NumWorkers > inferenceOptions.maxWorkers
                delete(gcp('nocreate')); % if pool with too many workers, delete and restart
                parpool(inferenceOptions.maxWorkers);
            end
            
            %% conduct cpHMM inference
            parfor i_local = 1:inferenceOptions.n_localEM % Parallel Local EM 
              
                % Random initialization of model parameters
                param_init = initialize_random_with_priors(inferenceOptions.nStates, noise_iid, v_iid);
                
                % Get Intial Values
                pi0_log_init = log(param_init.pi0);
                A_log_init = log(param_init.A);
                v_init = param_init.v;                        
                noise_init = param_init.noise;
                
                %--------------------LocalEM Call-------------------------%
                if ~inferenceOptions.truncInference
                  local_out = local_em_MS2_reduced_memory(fluo_data, ...
                      v_init, noise_init, pi0_log_init', A_log_init, inferenceOptions.nStates, inferenceOptions.nSteps, ...
                      alpha, inferenceOptions.nStinferenceOptions.epsMax, inferenceOptions.eps);                    
                else
                  local_out = local_em_MS2_reduced_memory_truncated(fluo_data, ...
                        v_init, noise_init, pi0_log_init', A_log_init, inferenceOptions.nStates, inferenceOptions.nSteps, ...
                    alpha, inferenceOptions.nStinferenceOptions.epsMax, inferenceOptions.eps);  
                end
                %---------------------------------------------------------%                
                % Save Results                 
                local_struct(i_local).subset_id = i_local;
                local_struct(i_local).logL = local_out.logL;                
                local_struct(i_local).A = exp(local_out.A_log);
                local_struct(i_local).v = exp(local_out.v_logs).*local_out.v_signs;
                local_struct(i_local).r = exp(local_out.v_logs).*local_out.v_signs / Tres;                                
                local_struct(i_local).noise = 1/exp(local_out.lambda_log);
                local_struct(i_local).pi0 = exp(local_out.pi0_log);
                local_struct(i_local).total_stinferenceOptions.eps = local_out.n_iter;               
                local_struct(i_local).soft_struct = local_out.soft_struct;               
            end
            
            %% Record output
            [~, max_index] = max([local_struct.logL]); % Get index of best result  
            
            % Save parameters from most likely local run
            output.pi0 =local_struct(max_index).pi0;                        
            output.r = local_struct(max_index).r(:);          
            output.noise = local_struct(max_index).noise;
            output.A = local_struct(max_index).A(:);
            output.A_mat = local_struct(max_index).A;  
            
            % get soft-decoded structure
            output.soft_struct = local_struct(max_index).soft_struct;
            
            % Info about run time
            output.total_stinferenceOptions.eps = local_struct(max_index).total_stinferenceOptions.eps;                                  
            output.total_time = 100000*(now - iter_start); 
            
            % other inference characteristics            
            output.protein_bin_flag = inferenceOptions.ProteinBinFlag;
            
            if inferenceOptions.ProteinBinFlag % NL: need to generalize this
                output.protein_bin = iter_index(t);
                output.protein_bin_list = iter_index;
                output.protein_bin_edges = mf_prctile_vec;
            end
                         
            output.iter_id = b;                        
            output.particle_ids = sample_particles;            
            output.N = ndp;
           
%             output.w = inferenceOptions.nSteps;
%             output.alpha = alpha;
%             output.deltaT = Tres;
%             output.inferenceOptions.fluo3DFlag = inferenceOptions.fluo3DFlag;
%             output.sampleSize = inferenceOptions.SampleSize; 
            % save inference data used
            output.fluo_data = fluo_data;
            output.time_data = time_data;
        end
        output.skip_flag = skip_flag;
        
        %% Determine unique filename and sace

        % Generate filenames            
        fName_sub = ['hmm_results_group' sprintf('%03d',t) '_rep'];
        file_list = dir([outDir fname_sub '*']);
        % Get largest sub-id
        if isempty(file_list) 
          repNum = 1;
        else
          repNumList = zeros(size(file_list));
          start_index = length(fname_sub)+1;
          for f = 1:length(file_list)
            repNumList(f) = str2double(file_list(f).name(start_index:start_index+2));
          end
          repNum = max(repNumList)+1;
        end
        % save
        out_file = [outDir '/' fName_sub sprintf('%03d',repNum)];          
        save([out_file '.mat'], 'output');           
    end  
end
 
