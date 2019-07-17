clear
close all

% define core ID variables
% project = 'Dl-Ven_snaBAC-mCh';
project = 'Dl-Ven_hbP2P-mCh';
dropboxFolder =  'E:\Nick\LivemRNA\Dropbox\';
dataPath = [dropboxFolder 'ProcessedEnrichmentData\' project '/'];
w = 7;
K = 3;
% window analysis params
window_size = 15; 
% arrays for fitting linear offsets
window_vec = -window_size:window_size;
fit_array = [ones(numel(window_vec),1) window_vec'];
% load input-output data set
load([dataPath 'hmm_input_output_w' num2str(w) '_K' num2str(K) '.mat'],'hmm_input_output')

% perform time averaging to get protein trends over time
gap_vec_full = [hmm_input_output.dt_filter_gap];
pt_vec_full = [hmm_input_output.spot_protein];
pt_vec_full = pt_vec_full(~gap_vec_full);
virtual_vec_full = [hmm_input_output.serial_protein];
virtual_vec_full = virtual_vec_full(~gap_vec_full);
time_vec_full_raw = [hmm_input_output.time];
time_vec_full = time_vec_full_raw(~gap_vec_full);

time_index = unique(time_vec_full_raw);
spot_protein_trend = NaN(size(time_index));
virtual_protein_trend = NaN(size(time_index));
t_sigma = 120;
for t = 1:numel(time_index)
    dt_vec = exp(-(time_vec_full-time_index(t)).^2 / 2 / t_sigma^2);
    % locus
    p_mean_spot = nansum(pt_vec_full.*dt_vec) ./ nansum(dt_vec);
    spot_protein_trend(t) = p_mean_spot;
    % virtual spot
    p_mean_virtual = nansum(virtual_vec_full.*dt_vec) ./ nansum(dt_vec);
    virtual_protein_trend(t) = p_mean_virtual;
end
% turn off display for fits
options = optimoptions('lsqnonlin','Display','off');
n_features = 0;
% iterate

for i = 1:numel(hmm_input_output)    
    z_vec = hmm_input_output(i).z_vec' > 1;
    r_vec = hmm_input_output(i).r_vec';
    hmm_input_output(i).z_vec_bin = z_vec;
    z_prob_vec = sum(hmm_input_output(i).z_mat(:,2:3),2);
    zd_full = [0 diff(z_vec)];
    change_points = find(zd_full~=0);
    rise_points = find(zd_full>0);
    % calculate duration of bursts
    burst_dist_vec_lead = diff([NaN rise_points]);
    dur_vec_lag = diff([change_points NaN]);
    dur_vec_lead = diff([NaN change_points]);
    % calculate intensity of bursts
    sz_vec_lag = NaN(size(dur_vec_lag));
    sz_vec_lead = NaN(size(dur_vec_lag));
    for j = 1:numel(change_points)-1
        r_mean = nanmean(r_vec(change_points(j):change_points(j+1)-1));
        sz_vec_lag(j) = r_mean;
        sz_vec_lead(j+1) = r_mean;
    end         
    % generate full-length vectors
    z_burst_dist_vec_full = NaN(size(z_vec));
    z_dur_lag_vec_full = NaN(size(z_vec));
    z_dur_lead_vec_full = NaN(size(z_vec));
    z_dur_lag_vec_full(change_points) = dur_vec_lag;  
    z_dur_lead_vec_full(change_points) = dur_vec_lead;  
    z_burst_dist_vec_full(rise_points) = burst_dist_vec_lead;
    sz_lag_vec_full = NaN(size(z_vec));
    sz_lead_vec_full = NaN(size(z_vec));
    sz_lag_vec_full(change_points) = sz_vec_lag;  
    sz_lead_vec_full(change_points) = sz_vec_lead; 

    % record
    hmm_input_output(i).z_burst_dist_vec = z_burst_dist_vec_full;  
    hmm_input_output(i).z_dur_lag_vec = z_dur_lag_vec_full;  
    hmm_input_output(i).z_dur_lead_vec = z_dur_lead_vec_full;
    hmm_input_output(i).sz_lag_vec = sz_lag_vec_full;  
    hmm_input_output(i).sz_lead_vec = sz_lead_vec_full;
    hmm_input_output(i).change_points = change_points;
    hmm_input_output(i).z_diff_vec = zd_full;
    hmm_input_output(i).z_prob_vec = z_prob_vec';
    % increment
    if ~isempty(change_points)
        n_features = n_features + numel(change_points)-1;
    end
    % generate detrended protein vectors
    dt_filter_gap = hmm_input_output(i).dt_filter_gap;
    dt_filter_gap_swap = hmm_input_output(i).dist_swap_dt_filter_gap;
    % locus protein       
    time_vec = hmm_input_output(i).time;
    % get corresponding pt trend points
    spot_trend_vec = spot_protein_trend(ismember(time_index,time_vec));
    virtual_trend_vec = virtual_protein_trend(ismember(time_index,time_vec));
    
    % perform simple fit to acount for differences in offset and dynamic
    % range
    spot_pt_vec = hmm_input_output(i).spot_protein; 
    spot_nan_ft = ~isnan(spot_pt_vec)&~dt_filter_gap;
    spot_protein_fit = detrend(spot_pt_vec(spot_nan_ft),2,'SamplePoints',time_vec(spot_nan_ft));
    spot_protein_dt = NaN(size(dt_filter_gap));
    spot_protein_dt(spot_nan_ft) = spot_protein_fit;
    hmm_input_output(i).spot_protein_dt = spot_protein_dt;

    % virtual spot protein    
    virtual_pt_vec = hmm_input_output(i).serial_protein;    
    nan_ft_virt = ~isnan(virtual_pt_vec)&~dt_filter_gap;
    serial_protein_fit = detrend(virtual_pt_vec(nan_ft_virt),2,'SamplePoints',time_vec(nan_ft_virt));     
    serial_protein_dt = NaN(size(dt_filter_gap));
    serial_protein_dt(nan_ft_virt) = serial_protein_fit;
    hmm_input_output(i).serial_protein_dt = serial_protein_dt;

    % swap spot protein
    swap_pt_vec = hmm_input_output(i).dist_swap_spot_protein;  
    nan_ft_swap = ~isnan(swap_pt_vec)&~dt_filter_gap_swap;
    swap_spot_protein_fit = detrend(swap_pt_vec(nan_ft_swap),2,'SamplePoints',time_vec(nan_ft_swap));
    swap_spot_protein_dt = NaN(size(dt_filter_gap));
    swap_spot_protein_dt(nan_ft_swap) = swap_spot_protein_fit;
    hmm_input_output(i).swap_spot_protein_dt = swap_spot_protein_dt;
end
%save update structure
save([dataPath 'hmm_input_output_w' num2str(w) '_K' num2str(K) '_dt.mat'],'hmm_input_output')
 
%%% Now compile snips for average burst dynamics analyses
% generate master set of vectors for feature classification
gap_filter_vec = [hmm_input_output.dt_filter_gap];
time_vec = round([hmm_input_output.time]/60);
z_dur_lag_list = [hmm_input_output.z_dur_lag_vec];
z_dur_lead_list = [hmm_input_output.z_dur_lead_vec];
sz_lag_list = [hmm_input_output.sz_lag_vec];
sz_lead_list = [hmm_input_output.sz_lead_vec];
z_diff_list = [hmm_input_output.z_diff_vec];
c = [hmm_input_output.z_diff_vec];
% initialize results structure
results_struct = struct;   
% calculate expected number of features for pre-allocation
n_entries = sum(z_diff_list~=0&~gap_filter_vec);
% initialize arrays to store time series snip
fluo_array = NaN(n_entries,2*window_size+1);
hmm_array = NaN(n_entries,2*window_size+1);
% de-meaned arrays
swap_hmm_array_dm = NaN(n_entries,2*window_size+1);
spot_array_dm = NaN(n_entries,2*window_size+1);
swap_array_dm = NaN(n_entries,2*window_size+1);
virtual_array_dm = NaN(n_entries,2*window_size+1);        
% detrended arrays
swap_hmm_array_dt = NaN(n_entries,2*window_size+1);
spot_array_dt = NaN(n_entries,2*window_size+1);
swap_array_dt = NaN(n_entries,2*window_size+1);
mf_array = NaN(n_entries,2*window_size+1);
virtual_array_dt = NaN(n_entries,2*window_size+1);        
% initialize arrays to store id variables
mf_protein_vec = NaN(1,n_entries);
particle_id_vec = NaN(1,n_entries);
center_time_vec = NaN(1,n_entries);
lag_dur_vec = NaN(1,n_entries);
lead_dur_vec = NaN(1,n_entries);
lag_size_vec = NaN(1,n_entries);
lead_size_vec = NaN(1,n_entries);
feature_sign_vec = NaN(1,n_entries);
prev_burst_dist_vec = NaN(1,n_entries);
% iterate through structure
iter = 1;                    
for j = 1:numel(hmm_input_output)
    % core ID vectors
    ParticleID = hmm_input_output(j).ParticleID;
    time = hmm_input_output(j).time;
    % feature classification vectors
    gap_filter = hmm_input_output(j).dt_filter_gap;   
    gap_filter_swap = hmm_input_output(j).dist_swap_dt_filter_gap;   
    z_dur_lag_vec_full = hmm_input_output(j).z_dur_lag_vec;            
    z_dur_lead_vec_full = hmm_input_output(j).z_dur_lead_vec;  
    sz_lag_vec_full = hmm_input_output(j).sz_lag_vec;            
    sz_lead_vec_full = hmm_input_output(j).sz_lead_vec;
    z_diff_vec = hmm_input_output(j).z_diff_vec;  
    z_burst_dist_vec = hmm_input_output(j).z_burst_dist_vec;
    % activity
    fluo = hmm_input_output(j).fluo;
    swap_r_vec = hmm_input_output(j).dist_swap_hmm;
    r_vec = hmm_input_output(j).r_vec';
    % raw protein fields
    spot_protein_raw = hmm_input_output(j).spot_protein;
    swap_spot_protein_raw = hmm_input_output(j).dist_swap_spot_protein;
    virtual_protein_raw = hmm_input_output(j).serial_protein;    
    mf_protein_raw = hmm_input_output(j).mf_protein;    
    % de-trended protein fields                        
    spot_protein_dt = hmm_input_output(j).spot_protein_dt;
    swap_spot_protein_dt = hmm_input_output(j).swap_spot_protein_dt;
    virtual_protein_dt = hmm_input_output(j).serial_protein_dt;        
    % apply filter to remove observations far too far from true points             
    spot_protein_dt(gap_filter) = NaN;
    spot_protein_raw(gap_filter) = NaN;
    swap_spot_protein_dt(gap_filter_swap) = NaN;
    swap_spot_protein_raw(gap_filter_swap) = NaN;
    virtual_protein_dt(gap_filter) = NaN;
    virtual_protein_raw(gap_filter) = NaN;
    mf_protein_raw(gap_filter) = NaN;
    % find features
    id_list = find(z_diff_vec~=0&~gap_filter);
    for id = id_list
        full_range = id - window_size:id+window_size;
        true_range = full_range(full_range>0&full_range<=numel(virtual_protein_dt));
        % record
        ft1 = ismember(full_range,true_range);
        ft2 = ~isnan(spot_protein_dt(true_range)) & ~isnan(swap_spot_protein_dt(true_range))&~isnan(virtual_protein_dt(true_range));
        % qc check
        if sum(ft2) >= window_size
            % extract rde-trended fragments 
            spot_fragment_dt = spot_protein_dt(true_range);% - nanmean(spot_protein(true_range));
            swap_fragment_dt = swap_spot_protein_dt(true_range);% - nanmean(swap_spot_protein(true_range));
            virtual_fragment_dt = virtual_protein_dt(true_range);% - nanmean(virtual_protein(true_range));
            % de-mean raw protein fragments
            spot_fragment_dm = spot_protein_raw(true_range) - nanmean(spot_protein_raw(true_range));
            swap_fragment_dm = swap_spot_protein_raw(true_range) - nanmean(swap_spot_protein_raw(true_range));
            virtual_fragment_dm = virtual_protein_raw(true_range) - nanmean(virtual_protein_raw(true_range));
            % output and mf fragments
            fluo_fragment = fluo(true_range);
            hmm_fragment = r_vec(true_range);
            swap_hmm_fragment = swap_r_vec(true_range);
            mf_fragment = mf_protein_raw(true_range);            
            % save de-trended protein snips               
            spot_array_dt(iter,ft1) = spot_fragment_dt;
            swap_array_dt(iter,ft1) = swap_fragment_dt;             
            virtual_array_dt(iter,ft1) = virtual_fragment_dt;
            % save de-meaned protein snips
            spot_array_dm(iter,ft1) = spot_fragment_dm;
            swap_array_dm(iter,ft1) = swap_fragment_dm;             
            virtual_array_dm(iter,ft1) = virtual_fragment_dm;
            % mf and output arrays
            mf_array(iter,ft1) = mf_fragment;
            fluo_array(iter,ft1) = fluo_fragment;
            hmm_array(iter,ft1) = hmm_fragment;
            swap_hmm_array_dt(iter,ft1) = swap_hmm_fragment;
            % save other info
            mf_protein_vec(iter) = nanmean(mf_fragment);
            particle_id_vec(iter) = ParticleID;
            center_time_vec(iter) = time(id);
            lag_dur_vec(iter) = z_dur_lag_vec_full(id);
            lead_dur_vec(iter) = z_dur_lead_vec_full(id);
            lag_size_vec(iter) = sz_lag_vec_full(id);
            lead_size_vec(iter) = sz_lead_vec_full(id);
            prev_burst_dist_vec(iter) = z_burst_dist_vec(id);
            feature_sign_vec(iter) = sign(z_diff_vec(id));
        end
        % increment
        iter = iter + 1;  
        if mod(iter,100) == 0
            disp(iter)
        end
    end    
end        
%%% record data
results_struct.mf_array = mf_array;
% detrended protein snips
results_struct.spot_array_dt = spot_array_dt;
results_struct.swap_array_dt = swap_array_dt;
results_struct.virtual_array_dt = virtual_array_dt;
% raw protein snips
results_struct.spot_array_dm = spot_array_dm;
results_struct.swap_array_dm = swap_array_dm;
results_struct.virtual_array_dm = virtual_array_dm;
% transcriptional response snips
results_struct.fluo_array = fluo_array;
results_struct.hmm_array = hmm_array;  
results_struct.swap_hmm_array = swap_hmm_array_dt;  
% ID variables
particle_id_vec(iter) = ParticleID;
results_struct.mf_protein_vec = mf_protein_vec;
results_struct.center_time_vec = center_time_vec;
results_struct.lag_dur_vec = lag_dur_vec;
results_struct.lead_dur_vec = lead_dur_vec;
results_struct.lag_size_vec = lag_size_vec;
results_struct.lead_size_vec = lead_size_vec;
results_struct.feature_sign_vec = feature_sign_vec;
results_struct.prev_burst_dist = prev_burst_dist_vec;
% save
save([dataPath 'hmm_input_output_results.mat'],'results_struct')