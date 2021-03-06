% Script to investigate distribution of enrichment "event" sizes 
clear
close all
addpath('utilities')
% define core ID variables
projectTarget = 'Dl-Ven_snaBAC-mCh';
projectControl = 'Dl-Ven_hbP2P-mCh';
dropboxFolder = 'E:\Nick\LivemRNA\Dropbox (Personal)\';
dataPath = [dropboxFolder 'ProcessedEnrichmentData\'];
figPath = [dropboxFolder 'LocalEnrichmentFigures\fluo_input_output\'];
mkdir(figPath)
% load data
w = 7;
K = 3;
load([dataPath projectTarget '\hmm_input_output_w' num2str(w) '_K' num2str(K) '.mat'],'hmm_input_output');
hmm_target = hmm_input_output;
load([dataPath projectControl '\hmm_input_output_w' num2str(w) '_K' num2str(K) '.mat'],'hmm_input_output');
hmm_control = hmm_input_output;
clear hmm_input_output
% Load analysis data
load([dataPath projectTarget '\nucleus_struct_protein.mat'], 'nucleus_struct_protein');
%% generate vectors of intensity changes and corresponding enrichment event size
% target
fluo_target = [hmm_target.fluo];
protein_delta_target = [hmm_target.spot_protein]-[hmm_target.null_protein];
mf_target = [hmm_target.mf_protein];

fluo_control = [hmm_control.fluo];
protein_delta_control = [hmm_control.spot_protein]-[hmm_control.null_protein];
mf_control = [hmm_control.mf_protein];

% analyze distributions
cm_hist = brewermap([],'Set2');

close all
pt_bins = linspace(-600,885,200);

overall_hist_pt = figure;
hold on
histogram(protein_delta_control,pt_bins,'Normalization','probability','FaceColor',cm_hist(3,:))
histogram(protein_delta_target,pt_bins,'Normalization','probability','FaceColor',cm_hist(2,:))
legend('enrichment (control)','enrichment (target)')
grid on
box on
set(gca,'Fontsize',14)
xlim([-500 500])
saveas(overall_hist_pt,[figPath 'full_enrichment_histograms.png'])
saveas(overall_hist_pt,[figPath 'full_enrichment_histograms.pdf'])


%% Fit distributions to enrichment data, filtered by spot fluorescence
close all
q_vec = [.2 .4 .6 .8];
fluo_q_target = quantile(fluo_target,q_vec);
% generate separate vectors for low, middle, and high spots
pt_target_dim = protein_delta_target(fluo_target<=fluo_q_target(1));
pt_target_mid = protein_delta_target(fluo_target>fluo_q_target(2)&fluo_target<=fluo_q_target(3));
pt_target_bright = protein_delta_target(fluo_target>fluo_q_target(4));

% calculate bootstrap averages and standard errors
n_boots = 100;
pt_target_array = NaN(n_boots,3);
for n = 1:n_boots
    pt_target_array(n,1) = nanmean(randsample(pt_target_dim,numel(pt_target_dim),true));
    pt_target_array(n,2) = nanmean(randsample(pt_target_mid,numel(pt_target_mid),true));
    pt_target_array(n,3) = nanmean(randsample(pt_target_bright,numel(pt_target_bright),true));
end
pt_target_mean = nanmean(pt_target_array);
pt_target_ste = nanstd(pt_target_array);

fluo_q_control = quantile(fluo_control,q_vec);
% generate separate vectors for low, middle, and high spots
pt_control_dim = protein_delta_control(fluo_control<=fluo_q_control(1));
pt_control_mid = protein_delta_control(fluo_control>fluo_q_control(2)&fluo_control<=fluo_q_control(3));
pt_control_bright = protein_delta_control(fluo_control>fluo_q_control(4));

pt_control_array = NaN(n_boots,3);
for n = 1:n_boots
    pt_control_array(n,1) = nanmean(randsample(pt_control_dim,numel(pt_control_dim),true));
    pt_control_array(n,2) = nanmean(randsample(pt_control_mid,numel(pt_control_mid),true));
    pt_control_array(n,3) = nanmean(randsample(pt_control_bright,numel(pt_control_bright),true));
end
pt_control_mean = nanmean(pt_control_array);
pt_control_ste = nanstd(pt_control_array);

%% Make filtered snip plots
DistLim = .8; % min distance from edge permitted (um)
PixelSize = nucleus_struct_protein(1).PixelSize;
dist_vec = [nucleus_struct_protein.spot_edge_dist_vec]*PixelSize;
dist_filter = true(size(dist_vec>=DistLim));

% fluorescence vector
fluo_vec = [nucleus_struct_protein.fluo];
% protein levels at center
spot_protein_vec = [nucleus_struct_protein.spot_protein_vec];
null_protein_vec = [nucleus_struct_protein.edge_null_protein_vec];
% Snip stacks
spot_protein_snips = cat(3,nucleus_struct_protein.spot_protein_snips);
null_protein_snips = cat(3,nucleus_struct_protein.edge_null_protein_snips);
spot_mcp_snips = cat(3,nucleus_struct_protein.spot_mcp_snips);
null_mcp_snips = cat(3,nucleus_struct_protein.edge_null_mcp_snips);

% apply distance filter (snips)
spot_protein_snips_dist = spot_protein_snips(:,:,dist_filter);
null_protein_snips_dist = null_protein_snips(:,:,dist_filter);
spot_mcp_snips_dist = spot_mcp_snips(:,:,dist_filter);
null_mcp_snips_dist = null_mcp_snips(:,:,dist_filter);
% (vectors)
fluo_vec_dist = fluo_vec(dist_filter);
spot_protein_vec_dist = spot_protein_vec(dist_filter); 
null_protein_vec_dist = null_protein_vec(dist_filter); 

% Make r reference array
snip_size = size(spot_protein_snips,1);
[y_ref, x_ref] = meshgrid(1:snip_size,1:snip_size);
r_ref = sqrt((x_ref-ceil(snip_size/2)).^2 + (y_ref-ceil(snip_size/2)).^2)*PixelSize;

inv_mat = [fliplr(1:snip_size); 1:snip_size]' ;
spot_protein_snips_mixed = NaN(size(spot_protein_snips_dist));
null_protein_snips_mixed = NaN(size(spot_protein_snips_dist));
spot_mcp_snips_mixed = NaN(size(spot_protein_snips_dist));
null_mcp_snips_mixed = NaN(size(spot_protein_snips_dist));
for i = 1:size(spot_protein_snips_dist,3)
    h = ceil(rand()*2);
    v = ceil(rand()*2);
    
    spot_protein_snips_mixed(:,:,i) = spot_protein_snips_dist(inv_mat(:,v),inv_mat(:,h),i);
    null_protein_snips_mixed(:,:,i) = null_protein_snips_dist(inv_mat(:,v),inv_mat(:,h),i);
    spot_mcp_snips_mixed(:,:,i) = spot_mcp_snips_dist(inv_mat(:,v),inv_mat(:,h),i);
    null_mcp_snips_mixed(:,:,i) = null_mcp_snips_dist(inv_mat(:,v),inv_mat(:,h),i);
end

%% Now  separate according to spot fluorescence
close all
% q_vec = [.2 .4 .6 .8];
fluo_q_target = quantile(fluo_vec_dist,q_vec);

% generate separate vectors for low, middle, and high spots
pt_target_dim = fluo_vec_dist<=fluo_q_target(1);
% pt_target_q2 = fluo_vec_dist>fluo_q_target(1)&fluo_vec_dist<=fluo_q_target(2);
pt_target_q3 = fluo_vec_dist>fluo_q_target(2)&fluo_vec_dist<=fluo_q_target(3);
pt_target_q4 = fluo_vec_dist>fluo_q_target(4);


% make average MCP maps
fluo_q1_snip_mean = nanmean(spot_mcp_snips_mixed(:,:,pt_target_dim),3);
fluo_q4_snip_mean = nanmean(spot_mcp_snips_mixed(:,:,pt_target_q4),3);

% make average protein snips
pt_q1_mean = nanmean(spot_protein_snips_mixed(:,:,pt_target_dim),3);
null_q1_mean = nanmean(null_protein_snips_mixed(:,:,pt_target_dim),3);
protein_q1_snip_mean = pt_q1_mean-null_q1_mean;

pt_q4_mean = nanmean(spot_protein_snips_mixed(:,:,pt_target_q4),3);
null_q4_mean = nanmean(null_protein_snips_mixed(:,:,pt_target_q4),3);
protein_q4_snip_mean = pt_q4_mean-null_q4_mean;

%% Save
px_size = nucleus_struct_protein(1).PixelSize;
z_size = nucleus_struct_protein(1).zStep;
voxel_size = z_size * px_size^2;

fluo_io_struct = struct;
fluo_io_struct.quantile_vec = q_vec;
% snips
fluo_io_struct.fluo_q1_snip = fluo_q1_snip_mean;
fluo_io_struct.fluo_q4_snip = fluo_q4_snip_mean;
fluo_io_struct.protein_q1_snip = protein_q1_snip_mean;
fluo_io_struct.protein_q4_snip = protein_q4_snip_mean;
% averages
fluo_io_struct.protein_target_mean = pt_target_mean;
fluo_io_struct.protein_target_ste = pt_target_ste;
fluo_io_struct.protein_control_mean = pt_control_mean;
fluo_io_struct.protein_control_ste = pt_control_ste;
% vox
fluo_io_struct.voxel_size = voxel_size;
% save
save([dataPath 'fluo_in_out.mat'],'fluo_io_struct')