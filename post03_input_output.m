% Script to probe relationship between input protein concentration and
% output transcriptional response
clear 
close all
% define ID variables
K = 3;
w = 7;
project = 'Dl-Ven x snaBAC';
nBoots = 100;
% project = 'Dl_Venus_hbP2P_MCPmCherry_Zoom2_7uW14uW';
dropboxFolder =  'E:\Nick\Dropbox (Garcia Lab)\';
% dropboxFolder = 'C:\Users\nlamm\Dropbox (Garcia Lab)\';
dataPath = [dropboxFolder '\ProcessedEnrichmentData\' project '\'];
figPath = [dropboxFolder '\LocalEnrichmentFigures\' project '\hmm_input_output_K' num2str(K) '_w' num2str(w) '\'];
mkdir(figPath)
load([dataPath 'input_output_snips.mat'])
gene_name = 'snaBAC';
protein_name = 'Dorsal';

%% Make time-dependent cross-covariance plots
% define some colors
yw = [234 194 100]/256; % yellow
bl = [115 143 193]/256; % blue
rd = [213 108 85]/256; % red
gr = [191 213 151]/256; % green
br = [207 178 147]/256; % brown

t_vec = [input_output_snips.t_center]; % center time for each snip
early_indices = find(t_vec < 20*60);
mid_indices = find(t_vec >= 60*20 & t_vec < 60*35);
late_indices = (t_vec >= 60*35);
n_lags = floor(numel(input_output_snips(3).ctrl_xcov)/2);

% conduct bootstrap sampling
early_xcov_mat = NaN(nBoots,2*n_lags+1);
mid_xcov_mat = NaN(nBoots,2*n_lags+1);
late_xcov_mat = NaN(nBoots,2*n_lags+1);
for n = 1:nBoots
    early_boot = randsample(early_indices,numel(early_indices),true);
    mid_boot = randsample(mid_indices,numel(mid_indices),true);
    late_boot = randsample(late_indices,numel(late_indices),true);

    early_xcov_mat(n,:) = nanmean(vertcat(input_output_snips(early_boot).spot_xcov));
    mid_xcov_mat(n,:) = nanmean(vertcat(input_output_snips(mid_boot).spot_xcov));
    late_xcov_mat(n,:) = nanmean(vertcat(input_output_snips(late_boot).spot_xcov));
end
early_xcov_mean = nanmean(early_xcov_mat);
early_xcov_ste = nanstd(early_xcov_mat);
mid_xcov_mean = nanmean(mid_xcov_mat);
mid_xcov_ste = nanstd(mid_xcov_mat);
late_xcov_mean = nanmean(late_xcov_mat);
late_xcov_ste = nanstd(late_xcov_mat);    

lag_axis = (-n_lags:n_lags)*20/60;
test_fig = figure;
hold on
fill([lag_axis fliplr(lag_axis)],[early_xcov_mean+early_xcov_ste fliplr(early_xcov_mean-early_xcov_ste)],bl,'FaceAlpha',.3,'EdgeAlpha',0)
fill([lag_axis fliplr(lag_axis)],[mid_xcov_mean+mid_xcov_ste fliplr(mid_xcov_mean-mid_xcov_ste)],gr,'FaceAlpha',.3,'EdgeAlpha',0)
fill([lag_axis fliplr(lag_axis)],[late_xcov_mean+late_xcov_ste fliplr(late_xcov_mean-late_xcov_ste)],rd,'FaceAlpha',.3,'EdgeAlpha',0)

p1 = plot(lag_axis,early_xcov_mean,'Color',bl);
p2 = plot(lag_axis,mid_xcov_mean,'Color',gr);
p3 = plot(lag_axis,late_xcov_mean,'Color',rd);
legend([p1 p2 p3],'early','middle','late')
grid on
xlabel('offset (minutes)')
ylabel(['cross-covariance (' gene_name ' x ' protein_name])
saveas(test_fig,[figPath 'temporal_xcov.png'])

%% now make feature plots
feature_cell = {'fluo_rise','fluo_fall','protein_peak','protein_trough'};
% initialize arrays
fluo_fluo_rise_mat = NaN(nBoots,2*n_lags+1);
fluo_fluo_fall_mat = NaN(nBoots,2*n_lags+1);
protein_fluo_rise_mat = NaN(nBoots,2*n_lags+1);
protein_fluo_fall_mat = NaN(nBoots,2*n_lags+1);
protein_protein_peak_mat = NaN(nBoots,2*n_lags+1);
fluo_protein_peak_mat = NaN(nBoots,2*n_lags+1);
protein_protein_trough_mat = NaN(nBoots,2*n_lags+1);
fluo_protein_trough_mat = NaN(nBoots,2*n_lags+1);
% get indices for sampling
protein_peak_indices = find([input_output_snips.pt_peak_flag]==1);% & [input_output_snips.pt_feature_prom] > 100);
protein_trough_indices = find([input_output_snips.pt_trough_flag]==1 & [input_output_snips.pt_feature_prom] > 100);
fluo_rise_indices = find([input_output_snips.fluo_change_flags]==1);% & [input_output_snips.fluo_feature_prom] > 50);
fluo_fall_indices = find([input_output_snips.fluo_change_flags]==-1);% & [input_output_snips.fluo_feature_prom] > 50);
%%
for n = 1:nBoots
    for i = 1:numel(feature_cell)
        % sample indices
        eval(['indices = ' feature_cell{i} '_indices;'])
        boot_indices = randsample(indices,numel(indices),true);
        % extract
        p_vec_spot = nanmean(vertcat(input_output_snips(boot_indices).spot_protein_vec));
        p_vec_swap = nanmean(vertcat(input_output_snips(boot_indices).swap_spot_protein_vec));
        f_vec_spot = nanmean(vertcat(input_output_snips(boot_indices).fluo_vec));
        f_vec_swap = nanmean(vertcat(input_output_snips(boot_indices).swap_fluo_vec));
        % record
        eval(['protein_' feature_cell{i} '_mat(n,:) = p_vec_spot - p_vec_swap;'])
        eval(['fluo_' feature_cell{i} '_mat(n,:) = f_vec_spot - f_vec_swap;'])
    end
end    
%%    
% & [input_output_snips.fluo_change_size] > 40;

for i = 1:numel(feature_cell)
    % record
    f_name = feature_cell{i};
    eval(['protein_' f_name '_mean = nanmean(protein_' f_name '_mat);'])
    eval(['protein_' f_name '_ste = nanstd(protein_' f_name '_mat);'])
    
    eval(['fluo_' f_name '_mean = nanmean(fluo_' f_name '_mat);'])
    eval(['fluo_' f_name '_ste = nanstd(fluo_' f_name '_mat);'])
end

dep_var_vec = {'p','p','f','f'};

for i = 1:numel(feature_cell)
    % record
    f_name = feature_cell{i};
    eval(['f_vec = fluo_' f_name '_mean;'])
    eval(['p_vec = protein_' f_name '_mean;'])
    fig = figure;
    hold on
    if strcmpi(dep_var_vec{i},'f')
        eval(['high = f_vec + fluo_' f_name '_ste;'])
        eval(['low = f_vec - fluo_' f_name '_ste;'])
        yyaxis left
        fill([lag_axis fliplr(lag_axis)],[high fliplr(low)],bl,'FaceAlpha',.2,'EdgeAlpha',0)
        p1 = plot(lag_axis,f_vec,'-','Color',bl,'LineWidth',1.3);
        yyaxis right
        p2 = plot(lag_axis,p_vec,'Color',rd,'LineWidth',1.3);
    else
        eval(['high = p_vec + protein_' f_name '_ste;'])
        eval(['low = p_vec - protein_' f_name '_ste;'])
        yyaxis right
        fill([lag_axis fliplr(lag_axis)],[high fliplr(low)],rd,'FaceAlpha',.2,'EdgeAlpha',0)
        p2 = plot(lag_axis,p_vec,'-','Color',rd,'LineWidth',1.3);
        yyaxis left
        p1 = plot(lag_axis,f_vec,'Color',bl,'LineWidth',1.3);
    end        
    legend([p1 p2],gene_name,protein_name)
    grid on
    xlabel('offset (minutes)')
end