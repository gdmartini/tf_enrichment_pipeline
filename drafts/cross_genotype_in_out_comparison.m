% Script to generate figures establishing presence of enrichment bursts at
% start of transcription bursts
clear 
% close all
addpath('../utilities')
% set ID variables
project1 = 'Dl-Ven_snaBAC-mCh_F-F-F_v1';
project2 = 'Dl-Ven_snaBAC-mCh_v3';
DropboxFolder = 'S:\Nick\Dropbox\';

% Params
fluo_dim = 3;
K = 3;
w = 7;

% set write paths
[~, DataPath1, FigureRoot] =   header_function(DropboxFolder, project1); 
[~, DataPath2, ~] =   header_function(DropboxFolder, project2); 

FigPath = [FigureRoot '\' project1 '\input_comparisons\'];
mkdir(FigPath)

% load data
% final results
load([DataPath1 'hmm_input_output_results_w' num2str(w) '_K' num2str(K) '_f' num2str(fluo_dim) '.mat'])
results_struct1 = results_struct;
load([DataPath2 'hmm_input_output_results_w' num2str(w) '_K' num2str(K) '_f' num2str(fluo_dim) '.mat'])
results_struct2 = results_struct;
clear results_struct;

% intermediate input/output set
load([DataPath1 'hmm_input_output_w' num2str(w) '_K' num2str(K) '_f' num2str(fluo_dim) 'D_dt.mat'],'hmm_input_output')
hmm1 = hmm_input_output;
load([DataPath2 'hmm_input_output_w' num2str(w) '_K' num2str(K) '_f' num2str(fluo_dim) 'D_dt.mat'],'hmm_input_output')
hmm2 = hmm_input_output;
clear nucleus_struct;

% raw compiled data
load([DataPath1 'nucleus_struct.mat'])
nc_struct1 = nucleus_struct;
load([DataPath2 'nucleus_struct.mat'])
nc_struct2 = nucleus_struct;
clear nucleus_struct;
%%
Tres = 20; % seconds
% extract relevant arrays from project 1
lag_dur_vec_1 = results_struct1.lag_dur_vec;
lead_dur_vec_1 = results_struct1.lead_dur_vec;
hmm_array_dm_1 = results_struct1.hmm_array;
hmm_array_dm_1 = hmm_array_dm_1 / nanstd(hmm_array_dm_1(:));
fluo_array_dm_1 = results_struct1.fluo_array  - nanmean(results_struct1.fluo_array,2);
fluo_array_dm_1 = fluo_array_dm_1 / nanstd(fluo_array_dm_1(:));
mf_array_dm_1 = results_struct1.mf_array - nanmean(results_struct1.mf_array,2);
mf_vec_1 = results_struct1.mf_protein_vec;
time_vec_1 = results_struct1.center_time_vec/60;
spot_array_dm_1 = results_struct1.spot_array_dm;
virtual_array_dm_1 = results_struct1.virtual_array_dm;
feature_sign_vec_1 = results_struct1.feature_sign_vec;
% now project 2
lag_dur_vec_2 = results_struct2.lag_dur_vec;
lead_dur_vec_2 = results_struct2.lead_dur_vec;
hmm_array_dm_2 = results_struct2.hmm_array;
hmm_array_dm_2 = hmm_array_dm_2 / nanstd(hmm_array_dm_2(:));
fluo_array_dm_2 = results_struct2.fluo_array  - nanmean(results_struct2.fluo_array,2);
fluo_array_dm_2 = fluo_array_dm_2 / nanstd(fluo_array_dm_2(:));
mf_array_dm_2 = results_struct2.mf_array - nanmean(results_struct2.mf_array,2);
mf_vec_2 = results_struct2.mf_protein_vec;
time_vec_2 = results_struct2.center_time_vec/60;
spot_array_dm_2 = results_struct2.spot_array_dm;
virtual_array_dm_2 = results_struct2.virtual_array_dm;
feature_sign_vec_2 = results_struct2.feature_sign_vec;

%  determine snip size
n_col = size(spot_array_dm_1,2);
window_size = floor(n_col/2);
time_axis = (-window_size:window_size)*Tres/60;

% set basic analyisis parameters
nBoots = 100; % number of bootstrap samples to use
% min_pause_len = 5; % minimum length of preceding OFF period (in time steps)
% max_pause_len = 10;
min_pause_len = 6; % minimum length of preceding OFF period (in time steps)
max_pause_len = 1000;
min_burst_len = 3;
max_burst_len = 1000;
% max_burst_len = 12;

% generate basic filter for target locus and computational controls
burst_ft_1 = feature_sign_vec_1 == 1&lead_dur_vec_1>=min_pause_len&lead_dur_vec_1<=max_pause_len...
    &lag_dur_vec_1>=min_burst_len&lag_dur_vec_1<=max_burst_len;%
burst_ft_2 = feature_sign_vec_2 == 1&lead_dur_vec_2>=min_pause_len&lead_dur_vec_2<=max_pause_len...
    &lag_dur_vec_2>=min_burst_len&lag_dur_vec_2<=max_burst_len;%
% sampling vectors
sample_options_1 = find(burst_ft_1);
sample_options_2 = find(burst_ft_2);


% (1) make de-trended input-output figure with controls
boot_spot_array_dm_1 = NaN(nBoots,n_col);
boot_virt_array_dm_1 = NaN(nBoots,n_col);
boot_hmm_array_dm_1 = NaN(nBoots,n_col);
boot_mf_array_dm_1 = NaN(nBoots,n_col);
boot_fluo_array_dm_1 = NaN(nBoots,n_col);
boot_spot_array_dm_2 = NaN(nBoots,n_col);
boot_virt_array_dm_2 = NaN(nBoots,n_col);
boot_hmm_array_dm_2 = NaN(nBoots,n_col);
boot_mf_array_dm_2 = NaN(nBoots,n_col);
boot_fluo_array_dm_2 = NaN(nBoots,n_col);
% take bootstrap samples
for n = 1:nBoots
    % project 1
    s_ids_1 = randsample(sample_options_1,numel(sample_options_1),true);        
    boot_spot_array_dm_1(n,:) = nanmean(spot_array_dm_1(s_ids_1,:));
    boot_virt_array_dm_1(n,:) = nanmean(virtual_array_dm_1(s_ids_1,:));
    boot_hmm_array_dm_1(n,:) = nanmean(hmm_array_dm_1(s_ids_1,:));
    boot_fluo_array_dm_1(n,:) = nanmean(fluo_array_dm_1(s_ids_1,:));
    boot_mf_array_dm_1(n,:) = nanmean(mf_array_dm_1(s_ids_1,:));
    % project 2
    s_ids_2 = randsample(sample_options_2,numel(sample_options_2),true);        
    boot_spot_array_dm_2(n,:) = nanmean(spot_array_dm_2(s_ids_2,:));
    boot_virt_array_dm_2(n,:) = nanmean(virtual_array_dm_2(s_ids_2,:));
    boot_hmm_array_dm_2(n,:) = nanmean(hmm_array_dm_2(s_ids_2,:));
    boot_fluo_array_dm_2(n,:) = nanmean(fluo_array_dm_2(s_ids_2,:));
    boot_mf_array_dm_2(n,:) = nanmean(mf_array_dm_2(s_ids_2,:));
end

% project #1
% calculate mean and standard error for spot
spot_mean_1 = nanmean(boot_spot_array_dm_1);
spot_ste_1 = nanstd(boot_spot_array_dm_1);
% calculate mean and standard error for virtual spot
virt_mean_1 = nanmean(boot_virt_array_dm_1);
virt_ste_1 = nanstd(boot_virt_array_dm_1);
% calculate mean and standard error for hmm trend
hmm_mean_1 = nanmean(boot_hmm_array_dm_1);
hmm_ste_1 = nanstd(boot_hmm_array_dm_1);
% calculate mean and standard error for raw fluoresecence
fluo_mean_1 = nanmean(boot_fluo_array_dm_1);
fluo_ste_1 = nanstd(boot_fluo_array_dm_1);
% calculate mean and standard error for mf protein
mf_mean_1 = nanmean(boot_mf_array_dm_1);
mf_ste_1 = nanstd(boot_mf_array_dm_1);

% project #2
% calculate mean and standard error for spot
spot_mean_2 = nanmean(boot_spot_array_dm_2);
spot_ste_2 = nanstd(boot_spot_array_dm_2);
% calculate mean and standard error for virtual spot
virt_mean_2 = nanmean(boot_virt_array_dm_2);
virt_ste_2 = nanstd(boot_virt_array_dm_2);
% calculate mean and standard error for hmm trend
hmm_mean_2 = nanmean(boot_hmm_array_dm_2);
hmm_ste_2 = nanstd(boot_hmm_array_dm_2);
% calculate mean and standard error for raw fluoresecence
fluo_mean_2 = nanmean(boot_fluo_array_dm_2);
fluo_ste_2 = nanstd(boot_fluo_array_dm_2);
% calculate mean and standard error for mf protein
mf_mean_2 = nanmean(boot_mf_array_dm_2);
mf_ste_2 = nanstd(boot_mf_array_dm_2);

%% make figure
cmap1 = brewermap([],'Set2');

burst_dt_comp_fig = figure;
% Dorsal activity
hold on
% locus (project 1)
p1 = plot(time_axis,spot_mean_1,'-','Color',cmap1(2,:),'LineWidth',2);
% locus (project 2)
p2 = plot(time_axis,spot_mean_2,'--','Color',cmap1(2,:),'LineWidth',2);
ax = gca;
ax.YColor = 'black';%cmap1(2,:);
% grid on
xlabel('offset (minutes)')
ylabel('relative Dl concentration (au)')
legend('Dl at {\it snail} locus (FFF)','Dl at {\it snail} locus (OG)', 'Location','northwest');
set(gca,'Fontsize',14,'xtick',-4:2:4)
chH = get(gca,'Children');
set(gca,'Children',flipud(chH));
set(gca,    'Box','off',...
            'Color',[228,221,209]/255,...            
            'TickLength',[0.02,0.05])    
burst_dt_comp_fig.Color = 'white';        
burst_dt_comp_fig.InvertHardcopy = 'off';
% save
saveas(burst_dt_comp_fig,[FigPath 'locus_trend_comparisons.tif'])
saveas(burst_dt_comp_fig,[FigPath 'locus_trend_comparisons.pdf'])

burst_dt_fig_1 = figure;
% Dorsal activity
hold on
% locus (project 1)
p1 = plot(time_axis,spot_mean_1,'-','Color',cmap1(2,:),'LineWidth',2);
% locus (project 2)
p2 = plot(time_axis,virt_mean_1,'--','Color',cmap1(3,:),'LineWidth',2);
ax = gca;
ax.YColor = 'black';%cmap1(2,:);
% grid on
xlabel('offset (minutes)')
ylabel('relative Dl concentration (au)')
legend('Dl at {\it snail} locus','Dl at control locus', 'Location','northwest');
set(gca,'Fontsize',14,'xtick',-4:2:4)
chH = get(gca,'Children');
set(gca,'Children',flipud(chH));
set(gca,    'Box','off',...
            'Color',[228,221,209]/255,...            
            'TickLength',[0.02,0.05])    
burst_dt_fig_1.Color = 'white';        
burst_dt_fig_1.InvertHardcopy = 'off';
% save
saveas(burst_dt_fig_1,[FigPath 'locus_trend_FFF.tif'])
saveas(burst_dt_fig_1,[FigPath 'locus_trend_FFF.pdf'])

burst_dt_fig_2 = figure;
% Dorsal activity
hold on
% locus (project 1)
p1 = plot(time_axis,spot_mean_2,'-','Color',cmap1(2,:),'LineWidth',2);
% locus (project 2)
p2 = plot(time_axis,virt_mean_2,'--','Color',cmap1(3,:),'LineWidth',2);
ax = gca;
ax.YColor = 'black';%cmap1(2,:);
% grid on
xlabel('offset (minutes)')
ylabel('relative Dl concentration (au)')
legend('Dl at {\it snail} locus','Dl at control locus', 'Location','northwest');
set(gca,'Fontsize',14,'xtick',-4:2:4)
chH = get(gca,'Children');
set(gca,'Children',flipud(chH));
set(gca,    'Box','off',...
            'Color',[228,221,209]/255,...            
            'TickLength',[0.02,0.05])    
burst_dt_fig_2.Color = 'white';        
burst_dt_fig_2.InvertHardcopy = 'off';
% save
saveas(burst_dt_fig_2,[FigPath 'locus_trend_OG.tif'])
saveas(burst_dt_fig_2,[FigPath 'locus_trend_OG.pdf'])

%% Check trends in other fields
hmm_comp_fig = figure;
% Dorsal activity
hold on
% locus (project 1)
p1 = plot(time_axis,hmm_mean_1,'-','Color',cmap1(3,:),'LineWidth',2);
% locus (project 2)
p2 = plot(time_axis,hmm_mean_2,'--','Color',cmap1(3,:),'LineWidth',2);
ax = gca;
ax.YColor = 'black';%cmap1(2,:);
% grid on
xlabel('offset (minutes)')
ylabel('hmm-decoded {\it snail} activity (au)')
legend('FFF','OG', 'Location','northwest');
set(gca,'Fontsize',14,'xtick',-4:2:4)
chH = get(gca,'Children');
set(gca,'Children',flipud(chH));
set(gca,    'Box','off',...
            'Color',[228,221,209]/255,...            
            'TickLength',[0.02,0.05])    
hmm_comp_fig.Color = 'white';        
hmm_comp_fig.InvertHardcopy = 'off';
% save
saveas(hmm_comp_fig,[FigPath 'locus_hmmtrend_comparisons.tif'])
saveas(hmm_comp_fig,[FigPath 'locus_hmm_trend_comparisons.pdf'])

fluo_comp_fig = figure;
% Dorsal activity
hold on
% locus (project 1)
p1 = plot(time_axis,fluo_mean_1,'-','Color',cmap1(4,:),'LineWidth',2);
% locus (project 2)
p2 = plot(time_axis,fluo_mean_2,'--','Color',cmap1(4,:),'LineWidth',2);
ax = gca;
ax.YColor = 'black';%cmap1(2,:);
% grid on
xlabel('offset (minutes)')
ylabel('raw {\it snail} fluorescence (au)')
legend('FFF','OG', 'Location','northwest');
set(gca,'Fontsize',14,'xtick',-4:2:4)
chH = get(gca,'Children');
set(gca,'Children',flipud(chH));
set(gca,    'Box','off',...
            'Color',[228,221,209]/255,...            
            'TickLength',[0.02,0.05])    
fluo_comp_fig.Color = 'white';        
fluo_comp_fig.InvertHardcopy = 'off';
% save
saveas(fluo_comp_fig,[FigPath 'locus_fluo_trend_comparisons.tif'])
saveas(fluo_comp_fig,[FigPath 'locus_fluo_trend_comparisons.pdf'])

mf_comp_fig = figure;
% Dorsal activity
hold on
% locus (project 1)
p1 = plot(time_axis,mf_mean_1,'-','Color',cmap1(5,:),'LineWidth',2);
% locus (project 2)
p2 = plot(time_axis,mf_mean_2,'--','Color',cmap1(5,:),'LineWidth',2);
ax = gca;
ax.YColor = 'black';%cmap1(2,:);
% grid on
xlabel('offset (minutes)')
ylabel('nuclear Dorsal (au)')
legend('FFF','OG', 'Location','northwest');
set(gca,'Fontsize',14,'xtick',-4:2:4)
chH = get(gca,'Children');
set(gca,'Children',flipud(chH));
set(gca,    'Box','off',...
            'Color',[228,221,209]/255,...            
            'TickLength',[0.02,0.05])    
mf_comp_fig.Color = 'white';        
mf_comp_fig.InvertHardcopy = 'off';
% save
saveas(mf_comp_fig,[FigPath 'locus_mf_trend_comparisons.tif'])
saveas(mf_comp_fig,[FigPath 'locus_mf_trend_comparisons.pdf'])

%% look at distribution of event times and protein
time_bins = 0:50;
time_fig = figure;
hold on
histogram(time_vec_1(burst_ft_1),time_bins,'Normalization','probability')
histogram(time_vec_2(burst_ft_2),time_bins,'Normalization','probability')
xlabel('event time (minutes)')
ylabel('share')
set(gca,'Fontsize',14)
legend('FFF','OG')
saveas(time_fig,[FigPath 'locus_event_time_comparisons.tif'])

dorsal_bins = linspace(.6,3);
mf_fig = figure;
hold on
histogram(mf_vec_1(burst_ft_1),dorsal_bins,'Normalization','probability')
histogram(mf_vec_2(burst_ft_2),dorsal_bins,'Normalization','probability')
xlabel('nuclear Dorsal (minutes)')
ylabel('share')
set(gca,'Fontsize',14)
legend('FFF','OG')
saveas(mf_fig,[FigPath 'locus_event_nuclear_protein_comparisons.tif'])

%% downsample OG to see if FFF noise can be explained by lack of statistics
boot_spot_array_ds_2 = NaN(nBoots,n_col);
boot_virt_array_ds_2 = NaN(nBoots,n_col);
rng(122);
% subset_ids_1 = randsample(sample_options_2,numel(sample_options_1),false);
% take bootstrap samples
for n = 1:nBoots
    % project 1
    s_ids_2 = randsample(sample_options_2,numel(sample_options_1),true);        
    boot_spot_array_ds_2(n,:) = nanmean(spot_array_dm_2(s_ids_2,:));
    boot_virt_array_ds_2(n,:) = nanmean(virtual_array_dm_2(s_ids_2,:));   
end


ds_spot_fig = figure;
% Dorsal activity
hold on
% locus (project 1)
p1 = plot(time_axis,boot_spot_array_ds_2','--','Color',[cmap1(2,:) .1],'LineWidth',1);
% locus (project 2)
p2 = plot(time_axis,spot_mean_1,'-','Color',cmap1(2,:),'LineWidth',2);
ax = gca;
ax.YColor = 'black';%cmap1(2,:);
% grid on
xlabel('offset (minutes)')
ylabel('relative Dl concentration (au)')
legend([p2 p1(1)],'Dl at {\it snail} locus (FFF)',...
    'Dl at {\it snail} locus (OG bootstraps)', 'Location','northeast');
set(gca,'Fontsize',14,'xtick',-4:2:4)
chH = get(gca,'Children');
set(gca,'Children',flipud(chH));
set(gca,    'Box','off',...
            'Color',[228,221,209]/255,...            
            'TickLength',[0.02,0.05])    
ds_spot_fig.Color = 'white';        
ds_spot_fig.InvertHardcopy = 'off';
% save
saveas(ds_spot_fig,[FigPath 'locus_trend_downsample_comparisons.tif'])
saveas(ds_spot_fig,[FigPath 'locus_trend_downsample_comparisons.pdf'])


ds_virt_fig = figure;
% Dorsal activity
hold on
% locus (project 1)
p1 = plot(time_axis,boot_virt_array_ds_2','--','Color',[cmap1(3,:) .1],'LineWidth',1);
% locus (project 2)
p2 = plot(time_axis,virt_mean_1,'-','Color',cmap1(3,:),'LineWidth',2);
ax = gca;
ax.YColor = 'black';%cmap1(2,:);
% grid on
xlabel('offset (minutes)')
ylabel('relative Dl concentration (au)')
legend([p2 p1(1)],'Dl at {\it snail} locus (FFF)',...
    'Dl at {\it snail} locus (OG bootstraps)', 'Location','northeast');
set(gca,'Fontsize',14,'xtick',-4:2:4)
chH = get(gca,'Children');
set(gca,'Children',flipud(chH));
set(gca,    'Box','off',...
            'Color',[228,221,209]/255,...            
            'TickLength',[0.02,0.05])    
ds_virt_fig.Color = 'white';        
ds_virt_fig.InvertHardcopy = 'off';
% save
saveas(ds_virt_fig,[FigPath 'virt_trend_downsample_comparisons.tif'])
saveas(ds_virt_fig,[FigPath 'virt_trend_downsample_comparisons.pdf'])

% spot_mean_ds_2 = nanmean(boot_spot_array_ds_2);
% virt_mean_ds_2 = nanmean(boot_virt_array_ds_2);
%% Get a sense for relative size of each set
n_nuclei = [numel(nc_struct2) numel(nc_struct1)];
n_frames = [nansum([nc_struct2.N]) nansum([nc_struct1.N])];
n_frames_qc = [nansum([nc_struct2([nc_struct2.qc_flag]==1).N]) nansum([nc_struct1([nc_struct1.qc_flag]==1).N])];
n_features = [nansum([hmm2.z_diff_vec]>0) nansum([hmm1.z_diff_vec]>0)]/5;
n_features_qc1 = [nansum([hmm2.z_diff_vec]>0&~[hmm2.dt_filter_gap]) nansum([hmm1.z_diff_vec]>0&~[hmm1.dt_filter_gap])]/5;
n_features_qc2 = [nansum(feature_sign_vec_2==1) nansum(feature_sign_vec_1==1)]/5;
n_features_final = [nansum(burst_ft_2) nansum(burst_ft_1)]/5;

% n_burst_features = [size([nc_struct2.qc_flag]) sum([nc_struct2.qc_flag])];
data_flow_array =  [n_frames ; n_frames_qc ; n_features ; n_features_qc1 ; n_features_qc2 ; n_features_final];

data_fig = figure;
plot(data_flow_array,'-o') 
set(gca,'YScale','Log')
grid on
ylim([100 1.05*max(data_flow_array(:))])
set(gca,'Fontsize',14)
ylabel('N datapoints')
xlabel('pipeline stage')
legend('OG','FFF')
saveas(data_fig,[FigPath 'data_flow_fig.tif'])

data_diff_fig = figure;
plot(data_flow_array./ data_flow_array(3,:),'-o') 
set(gca,'YScale','Log')
set(gca,'Fontsize',14)
ylabel('N datapoints (pegged)')
xlabel('pipeline stage')
legend('OG','FFF')
grid on
saveas(data_diff_fig,[FigPath 'data_flow_fig_pegged.tif'])
