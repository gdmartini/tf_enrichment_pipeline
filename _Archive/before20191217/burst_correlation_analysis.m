% Script to attempt a systematic dissection of various factors driving
% proteinxtranscription burst coincidence
clear
close all
addpath('utilities')
% define core ID variables
project = 'Dl-Ven_snaBAC-mCh';
% project = 'Dl-Ven_hbP2P-mCh';
dropboxFolder =  'E:\Nick\LivemRNA\Dropbox\';
dataPath = [dropboxFolder 'ProcessedEnrichmentData\' project '\'];
figPath = [dropboxFolder 'LocalEnrichmentFigures\' project '\burst_analyses\'];
mkdir(figPath)
% load data
load([dataPath 'hmm_input_output_results.mat'])
w = 7;
K = 3;
load([dataPath 'hmm_input_output_w' num2str(w) '_K' num2str(K) '.mat'],'hmm_input_output');

% define size of window of interest
roi_window = 6; 
window_size = 15;
start = window_size + 2;
% extract roi_vectors from wapo and locus arrays
locus_protein_vec = nanmean(results_struct.spot_array_dt(:,start:start + roi_window),2);
swap_protein_vec = nanmean(results_struct.swap_array_dt(:,start:start + roi_window),2);
mf_protein_vec = nanmean(results_struct.mf_array(:,start:start + roi_window),2);
% pull other trend vectors
feature_sign_vec = results_struct.feature_sign_vec';
lag_size_vec = results_struct.lag_size_vec';
lead_size_vec = results_struct.lead_size_vec';
lag_dur_vec = results_struct.lag_dur_vec';
lead_dur_vec = results_struct.lead_dur_vec';
tr_burst_size_vec = lag_dur_vec.*lag_size_vec;
% make rise filter
rise_ft = feature_sign_vec == 1;
analysis_ft = rise_ft & lead_dur_vec>5 & ~isnan(locus_protein_vec);
%%
%%% Make Bivariate heatmaps
% n_bins = 100;
close all
mf_bins = linspace(prctile(mf_protein_vec,5),prctile(mf_protein_vec,95),25);
burst_size_bins = linspace(prctile(tr_burst_size_vec,5),prctile(tr_burst_size_vec,95),25);
off_dur_bins = 1:11;
mf_sigma = median(diff(mf_bins))*2;
tr_size_sigma = median(diff(tr_burst_size_vec))*2;
dur_sigma = 1;
% amp_sigma = .1;
% initialize arrays
mf_curr_burst_mat = NaN(numel(mf_bins),numel(burst_size_bins));
mf_prev_burst_mat = NaN(numel(mf_bins),numel(off_dur_bins));
prev_curr_burst_mat = NaN(numel(burst_size_bins),numel(off_dur_bins));

% mean Dl vs. protein burst size
for m = 1:numel(mf_bins)    
    for d = 1:numel(burst_size_bins)
        dmf = (mf_bins(m) - mf_protein_vec(analysis_ft))/mf_sigma;
        ddu = (burst_size_bins(d) - tr_burst_size_vec(analysis_ft))/dur_sigma;
        wt = exp(-.5*(ddu.^2 + dmf.^2));
        mf_curr_burst_mat(m,d) = nansum(wt.*locus_protein_vec(analysis_ft)) / nansum(wt);
    end    
end
% make figure
mf_dur_fig = figure;
hm_cm = flipud(brewermap([],'RdYlBu'));
colormap(hm_cm);
imagesc(mf_curr_burst_mat);
h = colorbar;
ylabel(h,'protein burst size (au)')
ylabel('average Dl concentration')
xlabel('transcription burst size (au \times min)')
set(gca,'Fontsize',14)
set(gca,'ytick',1:4:25,'yticklabels',round(mf_bins(1:4:25)/5)*5)
set(gca,'xtick',1:4:25,'xticklabels',round(burst_size_bins(1:4:25)/3,1))
caxis([-10 40])
saveas(mf_dur_fig,[figPath 'mf_trSize_protein_hm.tif'])

% mean Dl vs. prev lag duration
for m = 1:numel(mf_bins)
    for d = 1:numel(off_dur_bins)
        dmf = (mf_bins(m) - mf_protein_vec(rise_ft))/mf_sigma;
        ddu = (off_dur_bins(d) - lead_dur_vec(rise_ft))/dur_sigma;
        wt = exp(-.5*(ddu.^2 + dmf.^2));
        mf_prev_burst_mat(m,d) = nansum(wt.*locus_protein_vec(rise_ft)) / nansum(wt);
    end
end

% make figure
prev_mf_fig = figure;
colormap(hm_cm);
imagesc(mf_prev_burst_mat');
h = colorbar;
ylabel(h,'protein burst size (au)')
xlabel('average Dl concentration (au)')
ylabel('preceding trough duration (min)')
set(gca,'Fontsize',14)
set(gca,'xtick',1:4:25,'xticklabels',round(mf_bins(1:4:25)/5)*5)
set(gca,'ytick',1:2:11,'yticklabels',round(off_dur_bins(1:2:11)/3,1))
caxis([-10 40])
saveas(prev_mf_fig,[figPath 'prev_mf_protein_hm.tif'])

% curr vs prev dur
for b = 1:numel(burst_size_bins)
    for d = 1:numel(off_dur_bins)
        db = (burst_size_bins(b) - tr_burst_size_vec(rise_ft))/dur_sigma;
        ddu = (off_dur_bins(d) - lead_dur_vec(rise_ft))/dur_sigma;
        wt = exp(-.5*(ddu.^2 + db.^2));
        prev_curr_burst_mat(b,d) = nansum(wt.*locus_protein_vec(rise_ft)) / nansum(wt);
    end
end

% make figure
prev_sz_fig = figure;
colormap(hm_cm);
imagesc(prev_curr_burst_mat');
h = colorbar;
ylabel(h,'protein burst size (au)')
xlabel('transcription burst size (au \times min)')
ylabel('preceding trough duration (min)')
set(gca,'Fontsize',14)
set(gca,'xtick',1:4:25,'xticklabels',round(burst_size_bins(1:4:25)/3,1))
set(gca,'ytick',1:2:11,'yticklabels',round(off_dur_bins(1:2:11)/3,1))
caxis([-10 40])
saveas(prev_sz_fig,[figPath 'prev_trSize_protein_hm.tif'])

%%
% % prev-mf heatmap
% amp_dur_fig = figure;
% colormap(hm_cm);
% p = imagesc(amp_dur_burst_mat);
% h = colorbar;
% ylabel(h,'protein burst size (au)')
% ylabel('transcription burst amplitude')
% xlabel('transcription burst duration')
% set(gca,'Fontsize',14)

%%% Examine protein burst duration as a function of transcription burst length
% spot_array = results_struct.spot_array;
% 
% spot_profile_mat = NaN(numel(burst_size_bins)-1, window_size+2);
% for i = 1:numel(burst_size_bins)-1
%     bdh = burst_size_bins(i+1);
%     if i+1 == numel(burst_size_bins)
%         bdh = Inf;
%     end
%     bdl = burst_size_bins(i);
%     iter_ft = rise_ft & lag_dur_vec >=bdl& lag_dur_vec <bdh & lead_dur_vec > 3;
%     spot_profile_mat(i,:) = nanmean(spot_array(iter_ft,window_size-2:end-2));
% end
% 
% % profile hm
% profile_fig = figure;
% colormap(hm_cm);
% p = imagesc(imgaussfilt(spot_profile_mat,1));
% h = colorbar;
% ylabel(h,'protein burst size (au)')
% ylabel('transcription burst duration')
% xlabel('position')
% set(gca,'Fontsize',14)
% % 
% % close all
% mf_index = 100:2:300;%linspace(prctile(mf_protein_vec(analysis_ft),5),prctile(mf_protein_vec(analysis_ft),95));
% pt_burst_index = -50:1.5:100;%linspace(prctile(locus_protein_vec(analysis_ft),5),prctile(locus_protein_vec(analysis_ft),99));
% 
% mf1 = mf_protein_vec(analysis_ft).^1;
% mf2 = mf_protein_vec(analysis_ft).^2;
% mf3 = mf_protein_vec(analysis_ft).^3;
% X = [mf1' mf2' mf3'];
% % fit to locus burst size
% mdl_locus = fitlm(X,locus_protein_vec(analysis_ft));
% mdl_swap = fitlm(X,swap_protein_vec(analysis_ft));
% 
% % generate predicted trends
% pd_burst_locus = predict(mdl_locus,[mf_index' (mf_index.^2)' (mf_index.^3)']);
% pd_burst_swap = predict(mdl_swap,[mf_index' (mf_index.^2)' (mf_index.^3)']);
% 
% % make figures
% locus_pd_scatter = figure;
% hold on
% scatter(mf_protein_vec(analysis_ft),locus_protein_vec(analysis_ft)',20,'MarkerFaceColor',...
%     [.6 .6 .6],'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',0);
% plot(mf_index,pd_burst_locus,'Color',hm_cm(end,:),'LineWidth',1.5)
% grid on
% xlabel('average Dl concentration')
% ylabel('protein burst size (locus)')
% legend('raw data','trend')
% set(gca,'Fontsize',14)
% xlim([100 300])
% ylim([-50 50])
% saveas(locus_pd_scatter,[figPath 'locus_mf_vs_pt_burst_scatter.tif'])
% 
% swap_pd_scatter = figure;
% hold on
% scatter(mf_protein_vec(analysis_ft),swap_protein_vec(analysis_ft)',20,'MarkerFaceColor',...
%     [.6 .6 .6],'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',0);
% plot(mf_index,pd_burst_swap,'Color',hm_cm(15,:),'LineWidth',1.5)
% grid on
% xlabel('average Dl concentration')
% ylabel('protein burst size (swap control)')
% legend('raw data','trend')
% set(gca,'Fontsize',14)
% xlim([100 300])
% ylim([-50 50])
% saveas(swap_pd_scatter,[figPath 'swap_mf_vs_pt_burst_scatter.tif'])
% 
% trend_comp_fig = figure;
% hold on
% plot(mf_index,pd_burst_locus,'Color',hm_cm(end,:),'LineWidth',1.5)
% plot(mf_index,pd_burst_swap,'Color',hm_cm(15,:),'LineWidth',1.5)
% grid on
% xlabel('average Dl concentration')
% ylabel('protein burst size (swap control)')
% legend('trend (locus)','trend (control)', 'Location','northwest')
% set(gca,'Fontsize',14)
% xlim([100 300])
% ylim([-5 20])
% saveas(trend_comp_fig,[figPath 'mf_vs_pt_burst_trends.tif'])










% %%
% rng(1); % For reproducibility
% Mdl = TreeBagger(100,[mf_protein_vec(analysis_ft)' locus_protein_vec(analysis_ft)],...
%     lag_dur_vec(analysis_ft)' .* lag_size_vec(analysis_ft)','Method','regression','OOBPredictorImportance','on');
% 
% yHat = oobPredict(Mdl);
% r2 = corr(Mdl.Y,yHat)^2;
% 
% MdlSwap = TreeBagger(100,[mf_protein_vec(analysis_ft)' swap_protein_vec(analysis_ft)],...
%     lag_dur_vec(analysis_ft)' .* lag_size_vec(analysis_ft)','Method','regression','OOBPredictorImportance','on');
% 
% yHatSwap = oobPredict(MdlSwap);
% r2Swap = corr(MdlSwap.Y,yHatSwap)^2;
% % mf_pt_scatter = figure;
% % colormap(hm_cm)
% % scatter(mf_protein_vec(analysis_ft),locus_protein_vec(analysis_ft),10,[.6 .6 .6],...
% %     'filled','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',0);
% 
% %%
% %%% Make Map of transcription burst duration as a function of local and average pt concentration
% close all
% local_burst_dur_mat = NaN(numel(mf_index),numel(pt_burst_index));
% swap_burst_dur_mat = NaN(numel(mf_index),numel(pt_burst_index));
% for m = 1:numel(mf_index)
%     for p = 1:numel(pt_burst_index)
%         mf = mf_index(m);
%         pt = pt_burst_index(p);
%         wt_vec_locus = exp(-.5*(((mf-mf_protein_vec(analysis_ft))/mf_sigma).^2+...
%             ((pt-locus_protein_vec(analysis_ft)')/pt_sigma).^2));
%         
%         wt_vec_swap = exp(-.5*(((mf-mf_protein_vec(analysis_ft))/mf_sigma).^2+...
%             ((pt-swap_protein_vec(analysis_ft)')/pt_sigma).^2));
%         
%         local_burst_dur_mat(m,p) = nansum(lag_dur_vec(analysis_ft).*wt_vec_locus) / nansum(wt_vec_locus);
%         swap_burst_dur_mat(m,p) = nansum(lag_dur_vec(analysis_ft).*wt_vec_swap) / nansum(wt_vec_swap);
%     end
% end
% 
% locus_mf_pt_fig = figure;
% colormap(hm_cm);
% p = imagesc(local_burst_dur_mat*20);
% h = colorbar;
% ylabel(h,'transcription burst duration (s)')
% ylabel('average Dl concentration')
% xlabel('local protein burst size')
% set(gca,'Fontsize',14)
% set(gca,'xtick',1:10:100,'xticklabels',round(pt_burst_index(1:10:100)))
% set(gca,'ytick',1:10:100,'yticklabels',round(mf_index(1:10:100)))
% caxis([2 9]*20)
% saveas(locus_mf_pt_fig,[figPath 'locus_mf_burst_dur_map.tif'])
% 
% 
% swap_mf_pt_fig = figure;
% colormap(hm_cm);
% p = imagesc(swap_burst_dur_mat*20);
% h = colorbar;
% ylabel(h,'transcription burst duration (s)')
% ylabel('average Dl concentration')
% xlabel('swap protein burst size')
% set(gca,'Fontsize',14)
% set(gca,'xtick',1:10:100,'xticklabels',round(pt_burst_index(1:10:100)))
% set(gca,'ytick',1:10:100,'yticklabels',round(mf_index(1:10:100)))
% caxis([2 9]*20)
% saveas(swap_mf_pt_fig,[figPath 'swap_mf_burst_dur_map.tif'])
% %%
% close all
% % prev burst duration
% prev_burst_scatter = figure;
% colormap(hm_cm)
% scatter(lead_dur_vec(rise_ft),locus_protein_vec(rise_ft),20,...
%     locus_protein_vec(rise_ft),'o','filled','MarkerFaceAlpha',.4,'MarkerEdgeAlpha',0)
% % ylim([0 30])
% xlim([0 30])
% grid on
% colorbar
% xlabel('preceding burst duration')
% ylabel('protein burst size')
% 
% 
% % coincident burst duration
% curr_burst_scatter = figure;
% colormap(hm_cm)
% scatter(lag_dur_vec(rise_ft),locus_protein_vec(rise_ft),20,...
%     locus_protein_vec(rise_ft),'o','filled','MarkerFaceAlpha',.4,'MarkerEdgeAlpha',0)
% % ylim([0 30])
% xlim([0 30])
% grid on
% colorbar
% xlabel('coincident burst duration')
% ylabel('protein burst size')
% 
% 
% % coincident burst duration
% mf_scatter = figure;
% colormap(hm_cm)
% scatter(mf_protein_vec(rise_ft),locus_protein_vec(rise_ft),20,...
%     locus_protein_vec(rise_ft),'o','filled','MarkerFaceAlpha',.4,'MarkerEdgeAlpha',0)
% % ylim([0 30])
% % xlim([0 30])
% grid on
% colorbar
% xlabel('average Dl concentration')
% ylabel('protein burst size')
% 
% 
% % coincident burst duration
% prev_amp_scatter = figure;
% colormap(hm_cm)
% scatter(lead_size_vec(rise_ft),locus_protein_vec(rise_ft),20,...
%     locus_protein_vec(rise_ft),'o','filled','MarkerFaceAlpha',.4,'MarkerEdgeAlpha',0)
% % ylim([0 30])
% % xlim([0 30])
% grid on
% colorbar
% xlabel('preceding burst amplitude')
% ylabel('protein burst size')
% 
% % coincident burst amplitude
% curr_amp_scatter = figure;
% colormap(hm_cm)
% scatter(lag_size_vec(rise_ft),locus_protein_vec(rise_ft),20,...
%     locus_protein_vec(rise_ft),'o','filled','MarkerFaceAlpha',.4,'MarkerEdgeAlpha',0)
% % ylim([0 30])
% % xlim([0 30])
% grid on
% colorbar
% xlabel('coincident burst amplitude')
% ylabel('protein burst size')
% 
% % coincident burst size
% curr_size_scatter = figure;
% colormap(hm_cm)
% scatter(hmm_curr_vec(rise_ft),locus_protein_vec(rise_ft),20,...
%     locus_protein_vec(rise_ft),'o','filled','MarkerFaceAlpha',.4,'MarkerEdgeAlpha',0)
% % ylim([0 30])
% % xlim([0 30])
% grid on
% colorbar
% xlabel('coincident burst size')
% ylabel('protein burst size')
% 
% 
% % prev burst size
% prev_size_scatter = figure;
% colormap(hm_cm)
% scatter(hmm_prev_vec(rise_ft),locus_protein_vec(rise_ft),20,...
%     locus_protein_vec(rise_ft),'o','filled','MarkerFaceAlpha',.4,'MarkerEdgeAlpha',0)
% % ylim([0 30])
% % xlim([0 30])
% grid on
% colorbar
% xlabel('preceding burst size')
% ylabel('protein burst size')
