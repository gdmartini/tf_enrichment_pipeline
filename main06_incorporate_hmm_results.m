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

function hmm_input_output = main06_incorporate_hmm_results(project,varargin)

close all
addpath('./utilities')
%%%%% These options will remain fixed for now
dpBootstrap = 0;
% dataRoot = ['../dat/'];
alphaFrac = 1302 / 6000;
dropboxFolder =  'E:\Nick\Dropbox (Garcia Lab)\';
dataPath = [dropboxFolder 'ProcessedEnrichmentData\' project '/'];
w = 7;
K = 3;

%%%%%%%%%%%%%%
for i = 1:numel(varargin)    
    if strcmpi(varargin{i},'dropboxFolder')
        dataRoot = [varargin{i+1} 'ProcessedEnrichmentData\'];
    end
    if ischar(varargin{i}) && i ~= numel(varargin)
        if ismember(varargin{i},{'dpBootstrap','controlProject'})
            eval([varargin{i} '=varargin{i+1};']);
        end
    end
end
d_type = '';
if dpBootstrap
    d_type = '_dp';
end
% extract 1c variables 
% dataPath = [dataRoot project '/'];
load([dataPath '/nucleus_struct_protein.mat'],'nucleus_struct_protein') % load data
load([dataPath '/nucleus_struct.mat'],'nucleus_struct') 
minDP = nucleus_struct_protein(1).minDP;
% check for necessary fields
analysis_fields = {'TresInterp','fluo_interp','time_interp'};
if true %~isfield(nucleus_struct_protein,analysis_fields{1})
    warning('Interpolation fields missing. Adding now.')    
    for i = 1:numel(nucleus_struct)
        for a = 1:numel(analysis_fields)
            nucleus_struct_protein(i).(analysis_fields{a}) = nucleus_struct(i).(analysis_fields{a});
        end
    end
end
Tres = nucleus_struct_protein(1).TresInterp; % Time Resolution
alpha = alphaFrac*w;
% generate alpha kernel 
alpha_kernel = NaN(1,w);
for i = 1:w
    if i < alpha
        alpha_kernel(i) = ((i-1) / alpha  + .5 * 1/alpha )*Tres;
    elseif i > alpha && (i-1) < alpha
        alpha_kernel(i) = Tres*(1 - .5*(alpha-i+1)*(1-(i-1)/alpha));
    else
        alpha_kernel(i) = Tres;
    end
end
voxel_size = nucleus_struct_protein(1).PixelSize^2 * nucleus_struct_protein(1).zStep;
% Set write path (inference results are now written to external directory)
hmm_suffix =  ['hmm_inference/w' num2str(w) '_K' num2str(K) '/']; 
file_list = dir([dataPath hmm_suffix 'hmm_results*.mat']);
if numel(file_list) > 1
    warning('multiple inference files detected. Ignoring all but first')
end

inference_results = load([dataPath hmm_suffix file_list(1).name]);
inference_results = inference_results.output;
% check for soft fit structure
soft_fit_flag = 1;
if exist([dataPath hmm_suffix 'soft_fit_struct.mat']) > 0
    fit_props = dir([dataPath hmm_suffix 'soft_fit_struct.mat']);
    fit_date = datenum(fit_props(1).date);
    hmm_date = datenum(file_list(1).date);
    if fit_date > hmm_date
        soft_fit_flag = 0;
    end
end
qc_indices = find([nucleus_struct_protein.qc_flag]==1);
particle_index = [nucleus_struct_protein.ParticleID];
if soft_fit_flag
    disp('conducting single trace fits...')
    A_log = log(inference_results.A_mat);
    v = inference_results.r*Tres;
    sigma = sqrt(inference_results.noise);
    pi0_log = log(inference_results.pi0); 
    eps = 1e-4;
    
    fluo_values = cell(numel(qc_indices),1);
    rm_indices = [];
    for i = 1:numel(qc_indices)
        fluo = nucleus_struct_protein(qc_indices(i)).fluo_interp;
        start_i = find(~isnan(fluo),1);
        stop_i = find(~isnan(fluo),1,'last');
        fluo = fluo(start_i:stop_i);
        if numel(fluo) < minDP
            error('problem with qc flag')
        end
        fluo_values{i} = fluo;
    end    
    fluo_values = fluo_values(~ismember(qc_indices,rm_indices));
    qc_indices = qc_indices(~ismember(qc_indices,rm_indices));
    tic 
    local_em_outputs = local_em_MS2_reduced_memory (fluo_values, ...
                  v', sigma, pi0_log, A_log, K, w, alpha, 1, eps);
    toc
    soft_fit_struct = local_em_outputs.soft_struct;
    soft_fit_struct.particle_index = particle_index(qc_indices);
    save([dataPath hmm_suffix 'soft_fit_struct.mat'],'soft_fit_struct')
else
    load([dataPath hmm_suffix 'soft_fit_struct.mat'],'soft_fit_struct')
end

%%% now extract corresponding hmm traces
hmm_input_output = [];
iter = 1;
for i = qc_indices%1:numel(particle_index)        
    % take average soft decoded result for particle
    ParticleID = particle_index(i);    
    if isnan(ParticleID)
        error('uh oh')
    end
    temp = struct;
    % extract relevant vectors from protein struct    
    ff_pt = nucleus_struct_protein(i).fluo;
    mcp_pt = nucleus_struct_protein(i).spot_mcp_vec;
    sp_pt = nucleus_struct_protein(i).spot_protein_vec;
    sp_pt_3D = nucleus_struct_protein(i).spot_protein_vec_3d;
    nn_pt = nucleus_struct_protein(i).edge_null_protein_vec;
    mf_pt = nucleus_struct_protein(i).mf_null_protein_vec;
    sr_pt = nucleus_struct_protein(i).serial_null_protein_vec;       
    sr_pt_3D = nucleus_struct_protein(i).serial_null_protein_vec_3d;
    tt_pt = nucleus_struct_protein(i).time;
    
    master_time = nucleus_struct_protein(i).time_interp;
    master_fluo = nucleus_struct_protein(i).fluo_interp;
    start_i = find(~isnan(master_fluo),1);
    stop_i = find(~isnan(master_fluo),1,'last');
    temp.time = master_time(start_i:stop_i);%inference_results(indices(1)).time_data{sub_indices(1)}; % doesn't matter which duplicate we reference
    temp.fluo = master_fluo(start_i:stop_i);%inference_results(indices(1)).fluo_data{sub_indices(1)};    
    % extract useful values
    [r,ri] = sort(inference_results.r);
    z = exp(soft_fit_struct.p_z_log_soft{iter});    
    temp.z_mat = z(ri,:)';    
    temp.r_mat = z(ri,:)'.*r';
    temp.r_inf = r';
    temp.r_vec = sum(temp.r_mat,2)';
    [~,z_vec] = max(temp.z_mat,[],2);
    temp.z_vec = z_vec; 
    % make predicted fluo vec
    fluo_hmm = conv(temp.r_vec,alpha_kernel);
    temp.fluo_hmm = fluo_hmm(1:numel(temp.r_vec));        
    if  sum(~isnan(sr_pt)&~isnan(mf_pt)) > 2
        % checks using mcp channel to ensure that we are correctly matching
        % particles and time frames
        temp.mcp_check = interp1(tt_pt(~isnan(mcp_pt)),mcp_pt(~isnan(mcp_pt)),temp.time);
        temp.fluo_check = interp1(tt_pt(~isnan(ff_pt)),ff_pt(~isnan(ff_pt)),temp.time);
        % record raw data vectors
        temp.spot_protein_raw = sp_pt;        
        temp.mf_protein_raw = mf_pt;
        temp.null_protein_raw = nn_pt;
        temp.serial_protein_raw = sr_pt;  
        temp.time_raw = tt_pt;  
        % protein information
        temp.spot_protein = interp1(tt_pt(~isnan(sp_pt)),sp_pt(~isnan(sp_pt)),temp.time);        
        temp.spot_protein_3D = interp1(tt_pt(~isnan(sp_pt_3D)),sp_pt_3D(~isnan(sp_pt_3D)),temp.time);        
        temp.mf_protein = interp1(tt_pt(~isnan(mf_pt)),mf_pt(~isnan(mf_pt)),temp.time);  
        temp.null_protein = interp1(tt_pt(~isnan(nn_pt)),nn_pt(~isnan(nn_pt)),temp.time);
        temp.serial_protein = interp1(tt_pt(~isnan(sr_pt)),sr_pt(~isnan(sr_pt)),temp.time);
        temp.serial_protein_3D = interp1(tt_pt(~isnan(sr_pt_3D)),sr_pt_3D(~isnan(sr_pt_3D)),temp.time);
        % reset values to NaN that are too far removed from true ref point        
        input_times = tt_pt(~isnan(ff_pt));
        dt_vec_gap = NaN(size(temp.time));
        for t = 1:numel(dt_vec_gap)
            dt_vec_gap(t) = min(abs(temp.time(t)-input_times));   
        end
        temp.dt_filter_gap = dt_vec_gap > 60;           
        % record general info for later use
        temp.ParticleID = ParticleID; 
        temp.Tres = Tres;
        hmm_input_output  = [hmm_input_output temp];
    end
    iter = iter + 1;
end
% find nearest neighbor particles
% generate array of average protein levels for each nucleus
time_vec = unique([nucleus_struct_protein.time_interp]);%nucleus_struct_protein(1).interpGrid;
mf_array = NaN(numel(time_vec),numel(hmm_input_output));
start_time_vec = NaN(size(hmm_input_output));
stop_time_vec = NaN(size(hmm_input_output));
set_vec = NaN(size(hmm_input_output));
for i = 1:numel(hmm_input_output)
    t_vec = hmm_input_output(i).time;
    mf_vec = hmm_input_output(i).mf_protein;
    mf_array(ismember(time_vec,t_vec),i) = mf_vec;
    start_time_vec(i) = min(t_vec);
    stop_time_vec(i) = max(t_vec);
    set_vec(i) = floor(hmm_input_output(i).ParticleID);
end
% now find closest match for each nucleus
for i = 1:numel(hmm_input_output)
    mf_i = mf_array(:,i);       
    dt_mf_vec = nanmean(abs(mf_array-mf_i));   
    setID = set_vec(i);
    % only accept traces that are as long or longer
    dt_mf_vec(start_time_vec>start_time_vec(i)|stop_time_vec<stop_time_vec(i)|set_vec~=setID) = NaN;
    dt_mf_vec(i) = NaN;
    [~, best_ind] = nanmin(dt_mf_vec);
    % record vales 
    time_swap = hmm_input_output(best_ind).time;  
    time_i = hmm_input_output(i).time; 
    % fill
    s_pt = NaN(size(time_i));
    s_pt(ismember(time_i,time_swap)) = hmm_input_output(best_ind).spot_protein(ismember(time_swap,time_i));
    sr_pt = NaN(size(time_i));
    sr_pt(ismember(time_i,time_swap)) = hmm_input_output(best_ind).serial_protein(ismember(time_swap,time_i));
    mf_pt = NaN(size(time_i));
    mf_pt(ismember(time_i,time_swap)) = hmm_input_output(best_ind).mf_protein(ismember(time_swap,time_i));
    s_fluo = NaN(size(time_i));
    s_fluo(ismember(time_i,time_swap)) = hmm_input_output(best_ind).fluo(ismember(time_swap,time_i));
    r_fluo = NaN(size(time_i));
    r_fluo(ismember(time_i,time_swap)) = hmm_input_output(best_ind).r_vec(ismember(time_swap,time_i));
    % record
    hmm_input_output(i).swap_ind = best_ind;
    hmm_input_output(i).swap_spot_protein = s_pt;
    hmm_input_output(i).swap_serial_protein = sr_pt;
    hmm_input_output(i).swap_mf_protein = mf_pt;
    hmm_input_output(i).swap_fluo = s_fluo;
    hmm_input_output(i).swap_hmm = r_fluo;
end

% save results
save([dataPath 'hmm_input_output_w' num2str(w) '_K' num2str(K) '.mat'],'hmm_input_output')