clear
close all

% define core ID variables
project = 'Dl-Ven_snaBAC-mCh';
% project = 'Dl-Ven_hbP2P-mCh';
dropboxFolder =  'E:\Nick\Dropbox (Garcia Lab)\';
dataPath = [dropboxFolder 'ProcessedEnrichmentData\' project '/'];
w = 7;
K = 3;

% window analysis params
window_size = 15; 
out_quantiles = 11;
% load input-output data set
load([dataPath 'hmm_input_output_w' num2str(w) '_K' num2str(K) '.mat'],'hmm_input_output')
% generate vectors indicating hmm changepoints
kernel_size = 1;
for i = 1:numel(hmm_input_output)    
    z_vec = hmm_input_output(i).z_vec' > 1;
    hmm_input_output(i).z_vec = z_vec;
    z_prob_vec = sum(hmm_input_output(i).z_mat(:,2:3),2);
    zd = [0 diff(z_vec)];
    change_points = find(zd~=0);
    dur_vec_lag = diff([change_points NaN]);
    dur_vec_lead = diff([NaN change_points]);
%     center_points = round((change_points(1:end-1) + change_points(2:end))/2 + rand(1,numel(change_points)-1)-.5);    
    z_dur_lag_vec = NaN(size(z_vec));
    z_dur_lead_vec = NaN(size(z_vec));
    z_dur_lag_vec(change_points) = dur_vec_lag;  
    z_dur_lead_vec(change_points) = dur_vec_lead;  
    hmm_input_output(i).z_dur_lag_vec = z_dur_lag_vec;  
    hmm_input_output(i).z_dur_lead_vec = z_dur_lead_vec;  
    hmm_input_output(i).z_diff_vec = zd;
    hmm_input_output(i).z_prob_vec = z_prob_vec';
end

gap_filter_vec = [hmm_input_output.dt_filter_gap];
time_vec = round([hmm_input_output.time]/60);
% z_ref = [hmm_input_output.z_vec];
% z_ref = 2*(z_ref-.5);
z_dur_lag_list = [hmm_input_output.z_dur_lag_vec];
z_dur_lead_list = [hmm_input_output.z_dur_lead_vec];
z_diff_list = [hmm_input_output.z_diff_vec];
% initialize results structure
results_struct = struct;
% establish time classes to investigate
time_class_cell = {1:60,1:15,16:25,26:35,36:45,46:60};
z_window = 2;
z_dur_sizes = 2:15;
z_dur_range_vec = [-fliplr(z_dur_sizes) z_dur_sizes];
min_sep_vec = [2 4 6];
% arrays for fitting linear offsets
window_vec = -window_size:window_size;
fit_array = [ones(numel(window_vec),1) window_vec'];
for ti = 1:numel(time_class_cell)    
    z_dur_key = [];
    z_sep_key = [];
    time_range = time_class_cell{ti};
    center_time = round(mean(time_range));
    time_filter = ismember(time_vec,time_range); 
    iter_filter = ~gap_filter_vec&time_filter;        
    % hmm 
    n_reps = numel(z_dur_range_vec)*numel(min_sep_vec);
    hmm_quantile_array = NaN(n_reps,2*window_size+1,out_quantiles);    
    % fluorescence
    fluo_quantile_array = NaN(n_reps,2*window_size+1,out_quantiles);    
    % target locus
    spot_quantile_array = NaN(n_reps,2*window_size+1,out_quantiles);    
    % swap locus
    swap_quantile_array = NaN(n_reps,2*window_size+1,out_quantiles);    
    % virtual spot
    virtual_quantile_array = NaN(n_reps,2*window_size+1,out_quantiles);        
    % iterate through feature list
    meta_iter = 1;
    for i = 1:numel(z_dur_range_vec)
        for k = 1:numel(min_sep_vec)
            feature_center = z_dur_range_vec(i);        
            feature_sign = sign(feature_center);
            if  feature_sign == 1
                key_list = z_dur_lag_list;
                sep_list = z_dur_lead_list;
            else
                sep_list = z_dur_lag_list;
                key_list = z_dur_lead_list;
            end
            z_values = abs(feature_center-z_window:feature_center+z_window);
            min_sep = min_sep_vec(k);
            % record lead and lag vales
            z_dur_key = [z_dur_key feature_center];
            z_sep_key = [z_sep_key min_sep];
    %         bd = z_dur_range_vec(i);
            n_entries = sum(ismember(key_list,z_values)&...
                sep_list>=min_sep&z_diff_list==feature_sign&~gap_filter_vec);
            if n_entries < 50
                continue
            end
            % initialize quantile arrays
            fluo_array = NaN(n_entries,2*window_size+1);
            hmm_array = NaN(n_entries,2*window_size+1);
            spot_array = NaN(n_entries,2*window_size+1);
            swap_array = NaN(n_entries,2*window_size+1);
            virtual_array = NaN(n_entries,2*window_size+1);        
            % iterate through structure
            iter = 1;                    
            for j = 1:numel(hmm_input_output)
                gap_filter = hmm_input_output(j).dt_filter_gap;
    %             z_state_vec = hmm_input_output(j).z_vec;
    %             z_state_vec = 2*(z_state_vec-.5);            
                z_dur_lag_vec = hmm_input_output(j).z_dur_lag_vec;            
                z_dur_lead_vec = hmm_input_output(j).z_dur_lead_vec;  
                diff_vec = hmm_input_output(j).z_diff_vec;
                if  feature_sign == 1
                    key_vec = z_dur_lag_vec;
                    sep_vec = z_dur_lead_vec;
                else
                    sep_vec = z_dur_lag_vec;
                    key_vec = z_dur_lead_vec;
                end
                fluo = hmm_input_output(j).fluo;
                z_vec = hmm_input_output(j).z_prob_vec';
                % protein fields                        
                spot_protein = hmm_input_output(j).spot_protein;
                swap_spot_protein = hmm_input_output(j).swap_spot_protein;
                virtual_protein = hmm_input_output(j).serial_protein;        
                % apply filter
                z_dur_lag_vec(gap_filter) = NaN;            
                spot_protein(gap_filter) = NaN;
                swap_spot_protein(gap_filter) = NaN;
                virtual_protein(gap_filter) = NaN;
                % find features
                id_list = find(ismember(key_vec,z_values)&sep_vec>=min_sep&feature_sign==diff_vec&~gap_filter);
                for id = id_list
                    full_range = id - window_size:id+window_size;
                    true_range = full_range(full_range>0&full_range<=numel(virtual_protein));
                    % record
                    ft1 = ismember(full_range,true_range);
                    if sum(~isnan(spot_protein)) >= window_size && sum(~isnan(swap_spot_protein)) >= window_size && sum(~isnan(virtual_protein)) >= window_size
                        spot_fragment = spot_protein(true_range);
                        swap_fragment = swap_spot_protein(true_range);
                        virtual_fragment = virtual_protein(true_range);
                        fluo_fragment = fluo(true_range);
                        hmm_fragment = z_vec(true_range);
                        % fit linear offsets            
                        fit_sub_array = fit_array(ft1,:);
                        spot_fit = fit_sub_array(~isnan(spot_fragment),:) \ spot_fragment(~isnan(spot_fragment))';
                        swap_fit = fit_sub_array(~isnan(swap_fragment),:) \ swap_fragment(~isnan(swap_fragment))';
                        virtual_fit = fit_sub_array(~isnan(virtual_fragment),:) \ virtual_fragment(~isnan(virtual_fragment))';                
                        fluo_fit = fit_sub_array(~isnan(fluo_fragment),:) \ fluo_fragment(~isnan(fluo_fragment))';
                        % save                
                        spot_array(iter,ft1) = spot_fragment - spot_fit(1) - spot_fit(2)*window_vec(ft1);
                        swap_array(iter,ft1) = swap_fragment - swap_fit(1) - swap_fit(2)*window_vec(ft1);
                        virtual_array(iter,ft1) = virtual_fragment - virtual_fit(1) - virtual_fit(2)*window_vec(ft1);
                        fluo_array(iter,ft1) = fluo_fragment - fluo_fit(1) - fluo_fit(2)*window_vec(ft1);
                        hmm_array(iter,ft1) = hmm_fragment;
                    end
                    % increment
                    iter = iter + 1;            
                end    
            end        
            % record quantiles
            fluo_quantile_array(meta_iter,:,:) = quantile(fluo_array,out_quantiles)';    
            hmm_quantile_array(meta_iter,:,:) = quantile(hmm_array,out_quantiles)';    
            spot_quantile_array(meta_iter,:,:) = quantile(spot_array,out_quantiles)';                   
            swap_quantile_array(meta_iter,:,:) = quantile(swap_array,out_quantiles)';               
            virtual_quantile_array(meta_iter,:,:) = quantile(virtual_array,out_quantiles)';    
            meta_iter = meta_iter + 1;
        end
    end    
    results_struct(ti).hmm_quantile_array = hmm_quantile_array;
    results_struct(ti).fluo_quantile_array = fluo_quantile_array;
    results_struct(ti).spot_quantile_array = spot_quantile_array;
    results_struct(ti).swap_quantile_array = swap_quantile_array;
    results_struct(ti).virtual_quantile_array = virtual_quantile_array;
    
    results_struct(ti).z_dur_key = z_dur_key;
    results_struct(ti).z_sep_key = z_sep_key;
    results_struct(ti).z_sep_key = z_sep_key;
%     results_struct(ti).hmm_mean_array = hmm_mean_array;
%     results_struct(ti).fluo_mean_array = fluo_mean_array;
%     results_struct(ti).spot_mean_array = spot_mean_array;
%     results_struct(ti).swap_mean_array = swap_mean_array;
%     results_struct(ti).virtual_mean_array = virtual_mean_array;
    results_struct(ti).time_range = time_range;
    results_struct(ti).window_vec = window_vec;
    results_struct(ti).z_dur_range_vec = z_dur_range_vec;    
    results_struct(ti).z_dur_size_vec = z_dur_sizes;    
end            

% save
save([dataPath 'hmm_input_output_results.mat'],'results_struct')