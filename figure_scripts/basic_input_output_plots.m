% Script to generate figures establishing presence of enrichment bursts at
% start of transcription bursts
clear 
close all
addpath('utilities')
% set ID variables
targetProject = 'Dl-Ven_snaBAC-mCh';
controlProject = 'Dl-Ven_hbP2P-mCh';
dropboxFolder =  'E:\Nick\LivemRNA\Dropbox\';
dataPathTarget = [dropboxFolder 'ProcessedEnrichmentData\' targetProject '\'];
dataPathControl = [dropboxFolder 'ProcessedEnrichmentData\' controlProject '\'];
figPath = [dropboxFolder 'LocalEnrichmentFigures\' targetProject '\_paper_figures\input_output01\'];
mkdir(figPath)

% load data
load([dataPathTarget 'hmm_input_output_results.mat'])
target_results_struct = results_struct;
load([dataPathControl 'hmm_input_output_results.mat'])
control_results_struct = results_struct;
Tres = 20; % seconds
% extract relevant arrays from target project 
lag_dur_vec_target = target_results_struct.lag_dur_vec;
lead_dur_vec_target = target_results_struct.lead_dur_vec;
hmm_array = target_results_struct.hmm_array;
swap_hmm_array = target_results_struct.swap_hmm_array;
spot_array_dt = target_results_struct.spot_array_dt;
spot_array_dm = target_results_struct.spot_array_dm;
swap_array_dt = target_results_struct.swap_array_dt;
virtual_array_dt = target_results_struct.virtual_array_dt;
feature_sign_vec_target = target_results_struct.feature_sign_vec;
% extract arrays from control project
lag_dur_vec_control = control_results_struct.lag_dur_vec;
lead_dur_vec_control = control_results_struct.lead_dur_vec;
biocontrol_array_dt = control_results_struct.spot_array_dt;
feature_sign_vec_control = control_results_struct.feature_sign_vec;
%  determine snip size
n_col = size(swap_array_dt,2);
window_size = floor(n_col/2);
time_axis = (-window_size:window_size)*Tres/60;

% set basic analyisis parameters
nBoots = 100; % number of bootstrap samples to use
min_pause_len = 5; % minimum length of preceding OFF period (in time steps)
min_burst_len = 2;

%%% (1) make basic input-output figure
close all
burst_rise_hmm_array = NaN(nBoots,n_col);
burst_rise_spot_array = NaN(nBoots,n_col);
% generate basic filter for target locus and computational controls
burst_ft_primary = feature_sign_vec_target == 1&lead_dur_vec_target>=min_pause_len&lag_dur_vec_target>min_burst_len; % filter for rise events
burst_ft_control = feature_sign_vec_control == 1&lead_dur_vec_control>=min_pause_len&lag_dur_vec_control>min_burst_len; % filter for rise events
sample_options_target = find(burst_ft_primary);
sample_options_control = find(burst_ft_control);
% take bootstrap samples
for n = 1:nBoots
    s_ids_target = randsample(sample_options_target,numel(sample_options_target),true);
    burst_rise_hmm_array(n,:) = nanmean(hmm_array(s_ids_target,:));
    burst_rise_spot_array(n,:) = nanmean(spot_array_dm(s_ids_target,:));
end
% calculate mean and standard error
burst_rise_hmm_mean = nanmean(burst_rise_hmm_array);
burst_rise_spot_mean = nanmean(burst_rise_spot_array);
burst_rise_spot_ste = nanstd(burst_rise_spot_array);
% calculate upper and lower bound vectors
br_spot_ub = burst_rise_spot_mean + burst_rise_spot_ste;
br_spot_lb = burst_rise_spot_mean - burst_rise_spot_ste;

% make figure
burst_trend_fig = figure;
cmap1 = brewermap([],'Set2');

% snail activity
yyaxis right
p1 = area(time_axis,burst_rise_hmm_mean,'FaceColor',[.6 .6 .6],'LineWidth',1.5,'FaceAlpha',.2);
ylabel('snail transcription (au)')
set(gca,'ytick',.2:.1:1.2)
ylim([.2 1.2])
ax = gca;
ax.YColor = 'black';

% Dorsal activity
yyaxis left
hold on
fill([time_axis fliplr(time_axis)],[br_spot_ub fliplr(br_spot_lb)],cmap1(2,:),'FaceAlpha',.2,'EdgeAlpha',0)
p2 = plot(time_axis,burst_rise_spot_mean,'-','Color',cmap1(2,:),'LineWidth',1.5);
ylabel('relative Dl concentration (au)')
set(gca,'ytick',-20:4:20)
ax = gca;
ax.YColor = 'black';

grid on
xlabel('offset (minutes)')
legend([p1 p2],'snail transcription','Dl concentration','Location','northwest')
set(gca,'Fontsize',12,'xtick',-4:2:4)
% save
saveas(burst_trend_fig,[figPath 'de-meaned_locus_trend.tif'])
saveas(burst_trend_fig,[figPath 'de-meaned_locus_trend.pdf'])

%%% (2) make de-trended input-output figure with controls

burst_rise_spot_array_dt = NaN(nBoots,n_col);
burst_rise_swap_array_dt = NaN(nBoots,n_col);
burst_rise_virt_array_dt = NaN(nBoots,n_col);
burst_rise_bio_array_dt = NaN(nBoots,n_col);
% take bootstrap samples
for n = 1:nBoots
    % primary
    s_ids_target = randsample(sample_options_target,numel(sample_options_target),true);    
    burst_rise_spot_array_dt(n,:) = nanmean(spot_array_dt(s_ids_target,:));
    burst_rise_swap_array_dt(n,:) = nanmean(swap_array_dt(s_ids_target,:));
    burst_rise_virt_array_dt(n,:) = nanmean(virtual_array_dt(s_ids_target,:));
    % biological control
    s_ids_control = randsample(sample_options_control,numel(sample_options_control),true);    
    burst_rise_bio_array_dt(n,:) = nanmean(biocontrol_array_dt(s_ids_control,:));
end

% calculate mean and standard error for spot
burst_rise_spot_mean = nanmean(burst_rise_spot_array_dt);
burst_rise_spot_ste = nanstd(burst_rise_spot_array_dt);
br_spot_ub = burst_rise_spot_mean + burst_rise_spot_ste;
br_spot_lb = burst_rise_spot_mean - burst_rise_spot_ste;
% calculate mean and standard error for nn swap
burst_rise_swap_mean = nanmean(burst_rise_swap_array_dt);
burst_rise_swap_ste = nanstd(burst_rise_swap_array_dt);
br_swap_ub = burst_rise_swap_mean + burst_rise_swap_ste;
br_swap_lb = burst_rise_swap_mean - burst_rise_swap_ste;
% calculate mean and standard error for virtual spot
burst_rise_virt_mean = nanmean(burst_rise_virt_array_dt);
burst_rise_virt_ste = nanstd(burst_rise_virt_array_dt);
br_virt_ub = burst_rise_virt_mean + burst_rise_virt_ste;
br_virt_lb = burst_rise_virt_mean - burst_rise_virt_ste;
% calculate mean and standard error for virtual spot
burst_rise_bio_mean = nanmean(burst_rise_bio_array_dt);
burst_rise_bio_ste = nanstd(burst_rise_bio_array_dt);
br_bio_ub = burst_rise_bio_mean + burst_rise_bio_ste;
br_bio_lb = burst_rise_bio_mean - burst_rise_bio_ste;

% make figure
burst_dt_fig = figure;
% snail activity
yyaxis right
p1 = area(time_axis,burst_rise_hmm_mean,'FaceColor',[.6 .6 .6],'LineWidth',1.5,'FaceAlpha',.2);
ylabel('snail transcription (au)')
set(gca,'ytick',.2:.1:1.2)
ylim([.2 1.2])
ax = gca;
ax.YColor = 'black';
% Dorsal activity
yyaxis left
hold on
% virtual control
fill([time_axis fliplr(time_axis)],[br_virt_ub fliplr(br_virt_lb)],cmap1(3,:),'FaceAlpha',.3,'EdgeAlpha',0)
p2 = plot(time_axis,burst_rise_virt_mean,'-','Color',cmap1(3,:),'LineWidth',2);
% swap control
fill([time_axis fliplr(time_axis)],[br_swap_ub fliplr(br_swap_lb)],cmap1(5,:),'FaceAlpha',.3,'EdgeAlpha',0)
p3 = plot(time_axis,burst_rise_swap_mean,'-','Color',cmap1(5,:),'LineWidth',2);
% biological control
fill([time_axis fliplr(time_axis)],[br_bio_ub fliplr(br_bio_lb)],cmap1(end,:),'FaceAlpha',.3,'EdgeAlpha',0)
p4 = plot(time_axis,burst_rise_bio_mean,'-','Color',cmap1(end,:),'LineWidth',2);
% locus
fill([time_axis fliplr(time_axis)],[br_spot_ub fliplr(br_spot_lb)],cmap1(2,:),'FaceAlpha',.5,'EdgeAlpha',0)
p5 = plot(time_axis,burst_rise_spot_mean,'-','Color',cmap1(2,:),'LineWidth',2);
ylabel('relative Dl concentration (au)')
set(gca,'ytick',-12:3:24)
ylim([-12 18])
ax = gca;
ax.YColor = 'black';
grid on
xlabel('offset (minutes)')
legend([p1 p2 p3 p4 p5],'snail transcription','virtual spot',...
    'nearest neighbor','off-target ({\ithbP2P})','target ({\itsnail})',...
    'Location','northwest')
set(gca,'Fontsize',12,'xtick',-4:2:4)
chH = get(gca,'Children');
set(gca,'Children',flipud(chH));
% save
saveas(burst_dt_fig,[figPath 'locus_trend_w_controls.tif'])
saveas(burst_dt_fig,[figPath 'locus_trend_w_controls.pdf'])

% make figure
burst_dt_fig_swap = figure;
% snail activity
yyaxis right
p1 = area(time_axis,burst_rise_hmm_mean,'FaceColor',[.6 .6 .6],'LineWidth',1.5,'FaceAlpha',.2);
ylabel('snail transcription (au)')
set(gca,'ytick',.2:.1:1.2)
ylim([.2 1.2])
ax = gca;
ax.YColor = 'black';
% Dorsal activity
yyaxis left
hold on
% swap control
fill([time_axis fliplr(time_axis)],[br_swap_ub fliplr(br_swap_lb)],cmap1(5,:),'FaceAlpha',.5,'EdgeAlpha',0)
p3 = plot(time_axis,burst_rise_swap_mean,'-','Color',cmap1(5,:),'LineWidth',2);
% locus
fill([time_axis fliplr(time_axis)],[br_spot_ub fliplr(br_spot_lb)],cmap1(2,:),'FaceAlpha',.5,'EdgeAlpha',0)
p5 = plot(time_axis,burst_rise_spot_mean,'-','Color',cmap1(2,:),'LineWidth',2);
ylabel('relative Dl concentration (au)')
set(gca,'ytick',-12:3:24)
ylim([-12 18])
ax = gca;
ax.YColor = 'black';
grid on
xlabel('offset (minutes)')
legend([p1 p3 p5],'snail transcription','control','target locus',...
    'Location','northwest')
set(gca,'Fontsize',12,'xtick',-4:2:4)
chH = get(gca,'Children');
set(gca,'Children',flipud(chH));
% save
saveas(burst_dt_fig,[figPath 'locus_trend_w_swap_control.tif'])
saveas(burst_dt_fig,[figPath 'locus_trend_w_swap_control.pdf'])

%%% (3) Make bar plots of preceding and succeeding 2 minutes
before_ids = time_axis < 0 & time_axis >=-2;
after_ids = time_axis > 0 & time_axis <=2;
% generate difference vectors
burst_rise_spot_after = nanmean(spot_array_dt(burst_ft_primary,after_ids),2);
burst_rise_spot_before = nanmean(spot_array_dt(burst_ft_primary,before_ids),2);

burst_rise_swap_after = nanmean(swap_array_dt(burst_ft_primary,after_ids),2);
burst_rise_swap_before = nanmean(swap_array_dt(burst_ft_primary,before_ids),2);

burst_rise_virt_after = nanmean(virtual_array_dt(burst_ft_primary,after_ids),2);
burst_rise_virt_before = nanmean(virtual_array_dt(burst_ft_primary,before_ids),2);

burst_rise_bio_after = nanmean(biocontrol_array_dt(burst_ft_control,after_ids),2);
burst_rise_bio_before = nanmean(biocontrol_array_dt(burst_ft_control,before_ids),2);

% Kolmogorov–Smirnov test)
[~, p_virt] = kstest2(burst_rise_virt_before,burst_rise_virt_after)
[~, p_swap] = kstest2(burst_rise_swap_before,burst_rise_swap_after)
[~, p_bio] = kstest2(burst_rise_bio_before,burst_rise_bio_after)
[~, p_spot] = kstest2(burst_rise_spot_before,burst_rise_spot_after)

bar_fig = figure;
cmap2 = brewermap([],'Paired');
colormap(cmap2);
bar_array = [nanmean(burst_rise_spot_before), nanmean(burst_rise_spot_after)
             nanmean(burst_rise_swap_before), nanmean(burst_rise_swap_after)
             nanmean(burst_rise_virt_before), nanmean(burst_rise_virt_after)
             nanmean(burst_rise_bio_before), nanmean(burst_rise_bio_after)];
         
b = bar(bar_array,'FaceColor','flat');
b(1).CData = cmap1([2 5 3 64],:);
b(1).FaceAlpha = .3;
b(2).CData = cmap1([2 5 3 64],:)/1.2;
set(gca,'xtick',1:4,'xticklabel',{'target locus','nearest neighbor','virtual spot','biological control',});
set(gca,'Fontsize',14)
xtickangle(30)
ylabel('Dorsal enrichment (au)')
ylim([-8 14])
grid on
saveas(bar_fig,[figPath 'bar_plots.tif'])
saveas(bar_fig,[figPath 'bar_plots.pdf'])