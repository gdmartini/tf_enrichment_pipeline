% Script to generate figures establishing presence of enrichment bursts at
% start of transcription bursts
clear 
close all
addpath('utilities')
% set ID variables
% targetProject = 'Dl-Ven_snaBAC-mCh_F-F-F_v1';
targetProject = '2xDl-Ven_snaBAC-mCh';
controlProject = '2xDl-Ven_hbP2P-mCh';
targetProjectOrig = 'Dl-Ven_snaBAC-mCh_v4';
project_cell = {targetProject controlProject};% targetProjectOrig};
DropboxFolder = 'S:\Nick\Dropbox\';
% DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\';
% Params
fluo_dim = 2;
protein_dim = 2;
K = 3;
w = 7;

[~, ~, FigureRoot] =   header_function(DropboxFolder, targetProject); 
FigPath = [FigureRoot '\' targetProject '\input_output01\'];
mkdir(FigPath)

% load data
master_struct = struct;
for i = 1:length(project_cell)
  % set write paths
  [~, DataPath] =   header_function(DropboxFolder, project_cell{i});   
  load([DataPath 'hmm_input_output_results_w' num2str(w) '_K' num2str(K) '_f' num2str(fluo_dim) 'D_p' num2str(protein_dim) 'D.mat'])
  master_struct(i).results_struct = results_struct;
  clear results_struct;
end


Tres = 20; % seconds
% extract relevant arrays from target project 
for i = 1:length(master_struct)
  master_struct(i).lag_dur_vec = master_struct(i).results_struct.lag_dur_vec;
  master_struct(i).lead_dur_vec = master_struct(i).results_struct.lead_dur_vec;
  master_struct(i).hmm_array = master_struct(i).results_struct.hmm_array;
  master_struct(i).swap_hmm_array = master_struct(i).results_struct.swap_hmm_array;
  master_struct(i).spot_array_dt = master_struct(i).results_struct.spot_array_dt;

  master_struct(i).swap_array_dt = master_struct(i).results_struct.swap_array_dt;
  master_struct(i).virtual_array_dt = master_struct(i).results_struct.virtual_array_dt;
  master_struct(i).feature_sign_vec = master_struct(i).results_struct.feature_sign_vec;
  
  % extract qc vectors 
  master_struct(i).swap_qc = master_struct(i).results_struct.swap_qc_vec;
  master_struct(i).virtual_qc = master_struct(i).results_struct.virtual_qc_vec;
  master_struct(i).set_vec = floor(master_struct(i).results_struct.particle_id_vec);
end

%  determine snip size
n_col = size(master_struct(1).swap_array_dt,2);
window_size = floor(n_col/2);
time_axis = (-window_size:window_size)*Tres/60;

% set basic analyisis parameters
nBoots = 100; % number of bootstrap samples to use
min_pause_len = 4; % minimum length of preceding OFF period (in time steps)
max_pause_len = 100;
min_burst_len = 2;
max_burst_len = 1000;


for i = 1:length(master_struct)
  % generate basic filter for target locus and computational controls
  master_struct(i).burst_ft = master_struct(i).feature_sign_vec == 1&master_struct(i).lead_dur_vec>=...
    min_pause_len&master_struct(i).lead_dur_vec<=max_pause_len...
      & master_struct(i).lag_dur_vec>=min_burst_len&master_struct(i).lag_dur_vec<=max_burst_len;%&target_swap_qc&target_virtual_qc;; % filter for rise events
  sample_options = find(master_struct(i).burst_ft);



  % (1) make de-trended input-output figure with controls
  master_struct(i).burst_rise_spot_array_dt = NaN(nBoots,n_col);
  master_struct(i).burst_rise_swap_array_dt = NaN(nBoots,n_col);
  master_struct(i).burst_rise_virt_array_dt = NaN(nBoots,n_col);
  master_struct(i).burst_rise_hmm_array = NaN(nBoots,n_col);
  % take bootstrap samples
  for n = 1:nBoots
      % primary
      s_ids_target = randsample(sample_options,numel(sample_options),true);    
      master_struct(i).burst_rise_hmm_array(n,:) = nanmean(master_struct(i).hmm_array(s_ids_target,:));
      master_struct(i).burst_rise_spot_array_dt(n,:) = nanmean(master_struct(i).spot_array_dt(s_ids_target,:));
      master_struct(i).burst_rise_swap_array_dt(n,:) = nanmean(master_struct(i).swap_array_dt(s_ids_target,:));
      master_struct(i).burst_rise_virt_array_dt(n,:) = nanmean(master_struct(i).virtual_array_dt(s_ids_target,:));      
  end

  % HMM trends
  master_struct(i).burst_rise_hmm_mean = nanmean(master_struct(i).burst_rise_hmm_array);
  master_struct(i).burst_rise_hmm_ste = nanstd(master_struct(i).burst_rise_hmm_array);
  % calculate mean and standard error for spot
  master_struct(i).burst_rise_spot_mean = nanmean(master_struct(i).burst_rise_spot_array_dt);
  master_struct(i).burst_rise_spot_ste = nanstd(master_struct(i).burst_rise_spot_array_dt);
  master_struct(i).br_spot_ub = master_struct(i).burst_rise_spot_mean + master_struct(i).burst_rise_spot_ste;
  master_struct(i).br_spot_lb = master_struct(i).burst_rise_spot_mean - master_struct(i).burst_rise_spot_ste;
  % calculate mean and standard error for nn swap
  master_struct(i).burst_rise_swap_mean = nanmean(master_struct(i).burst_rise_swap_array_dt);
  master_struct(i).burst_rise_swap_ste = nanstd(master_struct(i).burst_rise_swap_array_dt);
  master_struct(i).br_swap_ub = master_struct(i).burst_rise_swap_mean + master_struct(i).burst_rise_swap_ste;
  master_struct(i).br_swap_lb = master_struct(i).burst_rise_swap_mean - master_struct(i).burst_rise_swap_ste;
  % calculate mean and standard error for virtual spot
  master_struct(i).burst_rise_virt_mean = nanmean(master_struct(i).burst_rise_virt_array_dt);
  master_struct(i).burst_rise_virt_ste = nanstd(master_struct(i).burst_rise_virt_array_dt);
  master_struct(i).br_virt_ub = master_struct(i).burst_rise_virt_mean + master_struct(i).burst_rise_virt_ste;
  master_struct(i).br_virt_lb = master_struct(i).burst_rise_virt_mean - master_struct(i).burst_rise_virt_ste;
end
% calculate mean and standard error for virtual spot
% burst_rise_bio_mean = nanmean(burst_rise_bio_array_dt);
% burst_rise_bio_ste = nanstd(burst_rise_bio_array_dt);
% br_bio_ub = burst_rise_bio_mean + burst_rise_bio_ste;
% br_bio_lb = burst_rise_bio_mean - burst_rise_bio_ste;


%% make figure
cmap1 = brewermap([],'Set2');

burst_dt_fig_virt = figure;
% snail activity
yyaxis right
plot(time_axis,master_struct(1).burst_rise_hmm_mean,'--','LineWidth',2,'Color','black');
ylabel('snail transcription (au)')
ylim([.1 0.75])
% set(gca,'ytick',.1:.2:1.1)
ax = gca;
ax.YColor = 'black';

% Dorsal activity
yyaxis left
hold on

% % virtual control
% fill([time_axis fliplr(time_axis)],[master_struct(1).br_virt_ub fliplr(master_struct(1).br_virt_lb)],...
%   cmap1(3,:),'FaceAlpha',.2,'EdgeAlpha',0)
% p1 = plot(time_axis,master_struct(1).burst_rise_virt_mean,'-','Color',cmap1(3,:),'LineWidth',2);


% bio control
fill([time_axis fliplr(time_axis)],[master_struct(2).br_spot_ub fliplr(master_struct(2).br_spot_lb)],cmap1(3,:),'FaceAlpha',.2,'EdgeAlpha',0)
p3 = plot(time_axis,master_struct(2).burst_rise_spot_mean,'-','Color',cmap1(3,:),'LineWidth',2);
ylabel('relative Dl concentration (au)')

% locus
fill([time_axis fliplr(time_axis)],[master_struct(1).br_spot_ub fliplr(master_struct(1).br_spot_lb)],cmap1(2,:),'FaceAlpha',.5,'EdgeAlpha',0)
p2 = plot(time_axis,master_struct(1).burst_rise_spot_mean,'-','Color',cmap1(2,:),'LineWidth',2);
ylabel('relative Dl concentration (au)')

p = plot(0,0);

% ylim([-20 25])
% set(gca,'ytick',-20:5:25)
ax = gca;
ax.YColor = 'black';%cmap1(2,:);
% grid on
xlabel('offset (minutes)')
legend([p2 p3],'Dl at {\it snail} locus', 'Dl at {\it hbP2P}', 'Location','northwest');

set(gca,'Fontsize',14,'xtick',-4:2:4)
chH = get(gca,'Children');
set(gca,'Children',flipud(chH));
ylim([-.05 .07])
set(gca,    'Box','off',...
            'Color',[228,221,209]/255,...            
            'TickLength',[0.02,0.05])    
burst_dt_fig_virt.Color = 'white';        
burst_dt_fig_virt.InvertHardcopy = 'off';
% save
saveas(burst_dt_fig_virt,[FigPath 'snail_w_hbP2P_control.tif'])
saveas(burst_dt_fig_virt,[FigPath 'snail_w_hbP2P_control.pdf'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compare to original results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

burst_comparison = figure;
% snail activity
yyaxis right
plot(time_axis,master_struct(1).burst_rise_hmm_mean,'--','LineWidth',2,'Color','black');
ylabel('snail transcription (au)')
% ylim([.1 1.1])
% set(gca,'ytick',.1:.2:1.1)
ax = gca;
ax.YColor = 'black';

% Dorsal activity
yyaxis left
hold on

% % virtual control
% fill([time_axis fliplr(time_axis)],[master_struct(1).br_virt_ub fliplr(master_struct(1).br_virt_lb)],...
%   cmap1(3,:),'FaceAlpha',.5,'EdgeAlpha',0)
% plot(time_axis,master_struct(1).burst_rise_virt_mean,'-','Color',cmap1(3,:),'LineWidth',2);

% locus
fill([time_axis fliplr(time_axis)],[master_struct(1).br_spot_ub fliplr(master_struct(1).br_spot_lb)],cmap1(2,:),'FaceAlpha',.5,'EdgeAlpha',0)
p1 = plot(time_axis,master_struct(1).burst_rise_spot_mean,'-','Color',cmap1(2,:),'LineWidth',2);
ylabel('relative Dl concentration (au)')

% original result
fill([time_axis fliplr(time_axis)],[master_struct(3).br_spot_ub fliplr(master_struct(3).br_spot_lb)],brighten(cmap1(2,:),-0.5),'FaceAlpha',.5,'EdgeAlpha',0)
p2 = plot(time_axis,master_struct(3).burst_rise_spot_mean,'-','Color',brighten(cmap1(2,:),-0.5),'LineWidth',2);
ylabel('relative Dl concentration (au)')

p = plot(0,0);

% ylim([-20 25])
% set(gca,'ytick',-20:5:25)
ax = gca;
ax.YColor = 'black';%cmap1(2,:);
% grid on
xlabel('offset (minutes)')
legend([p1 p2],'Dl at {\it snail} (new)','Dl at {\it snail} (original)', 'Location','northwest');

set(gca,'Fontsize',14,'xtick',-4:2:4)
chH = get(gca,'Children');
set(gca,'Children',flipud(chH));
ylim([-.05 .1])
set(gca,    'Box','off',...
            'Color',[228,221,209]/255,...            
            'TickLength',[0.02,0.05])    
burst_comparison.Color = 'white';        
burst_comparison.InvertHardcopy = 'off';
% save
saveas(burst_comparison,[FigPath 'burst_comparison.tif'])
saveas(burst_comparison,[FigPath 'burst_comparison.pdf'])
