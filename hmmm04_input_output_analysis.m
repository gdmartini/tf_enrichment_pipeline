clear 
close all
% define ID variables
project = 'Dl_Venus_snaBAC_MCPmCherry_Zoom25x_minBleaching_test';
dropboxFolder =  'E:\Nick\Dropbox (Garcia Lab)\';
dataPath = [dropboxFolder '\ProcessedEnrichmentData\' project '\'];
K = 2;
w = 6;
n_lags = 10;
n_bins = 10;
% load data set
load([dataPath 'hmm_input_output_w' num2str(w) '_K' num2str(K) '.mat'])
% make burst vectors
rise_vec_full = [];
fall_vec_full = [];
r_vec_full = [];

particle_vec_full = [];

for i = 1:numel(hmm_input_output)
    rise_vec_full = [rise_vec_full [0 reshape(hmm_input_output(i).zz_mat(2,1,:),1,[])]];
    fall_vec_full = [fall_vec_full [0 reshape(hmm_input_output(i).zz_mat(1,2,:),1,[])]];
    r_vec_full = [r_vec_full reshape(sum(hmm_input_output(i).r_mat,2),1,[])];
    particle_vec_full = [particle_vec_full repelem(hmm_input_output(i).ParticleID,numel(hmm_input_output(i).time))];
end
time_vec_full = [hmm_input_output.time];
mf_vec = [hmm_input_output.mf_protein];
fluo_vec_full = [hmm_input_output.fluo];
delta_vec_full = [hmm_input_output.spot_protein] - [hmm_input_output.null_protein];

null_filter = isnan(delta_vec_full) | isnan(mf_vec);
time_vec_full = time_vec_full(~null_filter);
mf_vec = mf_vec(~null_filter);
delta_vec_full = delta_vec_full(~null_filter);
fluo_vec_full = fluo_vec_full(~null_filter);
rise_vec_full = rise_vec_full(~null_filter);
fall_vec_full = fall_vec_full(~null_filter);
r_vec_full = r_vec_full(~null_filter);
particle_vec_full = particle_vec_full(~null_filter);

set_vec_full = floor(particle_vec_full);

% identify outliers
low_ids = find(r_vec_full<prctile(r_vec_full,20));
high_ids = find(r_vec_full>prctile(r_vec_full,80));
rise_ids = find(rise_vec_full>.6);
fall_ids = find(fall_vec_full>.6);

% calculate sampling weights that allow us to draw control distributions
% that mimic mf and fluo values for subsets of interest
index_vec = 1:numel(fluo_vec_full);
% generate ID vector assigning each observation to a bin in 2D array
% bin_id_vec = NaN(size(fluo_vec));

% for i = 1:n_bins
%     for j = 1:n_bins
%         ft = fluo_vec < Xedges(i+1) & fluo_vec >= Xedges(i) & mf_vec < Yedges(j+1) & mf_vec >= Yedges(j);
%         bin_id_vec(ft) = sub2ind([n_bins,n_bins],j,i);
%     end
% end
% now calculate resampling weights fior each scenario
low_ft = ismember(index_vec,low_ids);
[baseHist,Xedges,Yedges,binX,binY] = histcounts2(fluo_vec_full(~low_ft),mf_vec(~low_ft),n_bins,'Normalization','probability');
bin_id_vec = sub2ind([n_bins,n_bins],binX,binY);
baseHist = baseHist + 1e-6;

lowHist = histcounts2(fluo_vec_full(low_ids),mf_vec(low_ids),Xedges,Yedges,'Normalization','probability');
low_wt = lowHist ./ baseHist;
low_wt_vec = low_wt(bin_id_vec);
low_ctrl_ids = randsample(index_vec(~low_ft),numel(low_ids),true,low_wt_vec);

% high
high_ft = ismember(index_vec,high_ids);
[baseHist,Xedges,Yedges,binX,binY] = histcounts2(fluo_vec_full(~high_ft),mf_vec(~high_ft),n_bins,'Normalization','probability');
bin_id_vec = sub2ind([n_bins,n_bins],binX,binY);
baseHist = baseHist + 1e-6;

highHist = histcounts2(fluo_vec_full(high_ids),mf_vec(high_ids),Xedges,Yedges,'Normalization','probability');
high_wt = highHist ./ baseHist;
high_wt_vec = high_wt(bin_id_vec);
high_ctrl_ids = randsample(index_vec(~high_ft),numel(high_ids),true,high_wt_vec);

% rise
rise_ft = ismember(index_vec,rise_ids);
[baseHist,Xedges,Yedges,binX,binY] = histcounts2(fluo_vec_full(~rise_ft),mf_vec(~rise_ft),n_bins,'Normalization','probability');
bin_id_vec = sub2ind([n_bins,n_bins],binX,binY);
baseHist = baseHist + 1e-6;

riseHist = histcounts2(fluo_vec_full(rise_ids),mf_vec(rise_ids),Xedges,Yedges,'Normalization','probability');
rise_wt = riseHist ./ baseHist;
rise_wt_vec = rise_wt(bin_id_vec);
rise_ctrl_ids = randsample(index_vec(~rise_ft),numel(rise_ids),true,rise_wt_vec);


% fall 
fall_ft = ismember(index_vec,fall_ids);
[baseHist,Xedges,Yedges,binX,binY] = histcounts2(fluo_vec_full(~fall_ft),mf_vec(~fall_ft),n_bins,'Normalization','probability');
bin_id_vec = sub2ind([n_bins,n_bins],binX,binY);
baseHist = baseHist + 1e-6;

fallHist = histcounts2(fluo_vec_full(fall_ids),mf_vec(fall_ids),Xedges,Yedges,'Normalization','probability');
fall_wt = fallHist ./ baseHist;
fall_wt_vec = fall_wt(bin_id_vec);
fall_ctrl_ids = randsample(index_vec(~fall_ft),numel(fall_ids),true,fall_wt_vec);

% low
low_array = NaN(numel(low_ids),n_lags+1);
low_ctrl_array = NaN(numel(low_ids),n_lags+1);
for i = 1:numel(low_ids)
    % sample
    ParticleID = particle_vec_full(low_ids(i));
    pt_vec = [NaN(1,max(0,n_lags-low_ids(i)+1)) particle_vec_full(max(1,low_ids(i)-n_lags):low_ids(i))];
    dt_vec = [NaN(1,max(0,n_lags-low_ctrl_ids(i)+1)) delta_vec_full(max(1,low_ctrl_ids(i)-n_lags):low_ctrl_ids(i))];
    low_array(i,pt_vec==ParticleID) = dt_vec(pt_vec==ParticleID);
    % control
    ParticleID = particle_vec_full(low_ctrl_ids(i));
    pt_vec = [NaN(1,max(0,n_lags-low_ctrl_ids(i)+1)) particle_vec_full(max(1,low_ctrl_ids(i)-n_lags):low_ctrl_ids(i))];
    dt_vec = [NaN(1,max(0,n_lags-low_ctrl_ids(i)+1)) delta_vec_full(max(1,low_ctrl_ids(i)-n_lags):low_ctrl_ids(i))];
    low_ctrl_array(i,pt_vec==ParticleID) = dt_vec(pt_vec==ParticleID);
end


high_array = NaN(numel(high_ids),n_lags+1);
high_ctrl_array = NaN(numel(high_ids),n_lags+1);
for i = 1:numel(high_ids)
    % sample
    ParticleID = particle_vec_full(high_ids(i));
    pt_vec = [NaN(1,max(0,n_lags-high_ids(i)+1)) particle_vec_full(max(1,high_ids(i)-n_lags):high_ids(i))];
    dt_vec = [NaN(1,max(0,n_lags-high_ctrl_ids(i)+1)) delta_vec_full(max(1,high_ctrl_ids(i)-n_lags):high_ctrl_ids(i))];
    high_array(i,pt_vec==ParticleID) = dt_vec(pt_vec==ParticleID);
    % control
    ParticleID = particle_vec_full(high_ctrl_ids(i));
    pt_vec = [NaN(1,max(0,n_lags-high_ctrl_ids(i)+1)) particle_vec_full(max(1,high_ctrl_ids(i)-n_lags):high_ctrl_ids(i))];
    dt_vec = [NaN(1,max(0,n_lags-high_ctrl_ids(i)+1)) delta_vec_full(max(1,high_ctrl_ids(i)-n_lags):high_ctrl_ids(i))];
    high_ctrl_array(i,pt_vec==ParticleID) = dt_vec(pt_vec==ParticleID);
end


rise_array = NaN(numel(rise_ids),n_lags+1);
rise_ctrl_array = NaN(numel(rise_ids),n_lags+1);

for i = 1:numel(rise_ids)
    % sample
    ParticleID = particle_vec_full(rise_ids(i));
    pt_vec = [NaN(1,max(0,n_lags-rise_ids(i)+1)) particle_vec_full(max(1,rise_ids(i)-n_lags):rise_ids(i))];
    dt_vec = [NaN(1,max(0,n_lags-rise_ctrl_ids(i)+1)) delta_vec_full(max(1,rise_ctrl_ids(i)-n_lags):rise_ctrl_ids(i))];
    rise_array(i,pt_vec==ParticleID) = dt_vec(pt_vec==ParticleID);
    % control
    ParticleID = particle_vec_full(rise_ctrl_ids(i));
    pt_vec = [NaN(1,max(0,n_lags-rise_ctrl_ids(i)+1)) particle_vec_full(max(1,rise_ctrl_ids(i)-n_lags):rise_ctrl_ids(i))];
    dt_vec = [NaN(1,max(0,n_lags-rise_ctrl_ids(i)+1)) delta_vec_full(max(1,rise_ctrl_ids(i)-n_lags):rise_ctrl_ids(i))];
    rise_ctrl_array(i,pt_vec==ParticleID) = dt_vec(pt_vec==ParticleID);
end


fall_array = NaN(numel(fall_ids),n_lags+1);
fall_ctrl_array = NaN(numel(fall_ids),n_lags+1);

for i = 1:numel(fall_ids)
    % sample
    ParticleID = particle_vec_full(fall_ids(i));
    pt_vec = [NaN(1,max(0,n_lags-fall_ids(i)+1)) particle_vec_full(max(1,fall_ids(i)-n_lags):fall_ids(i))];
    dt_vec = [NaN(1,max(0,n_lags-fall_ctrl_ids(i)+1)) delta_vec_full(max(1,fall_ctrl_ids(i)-n_lags):fall_ctrl_ids(i))];
    fall_array(i,pt_vec==ParticleID) = dt_vec(pt_vec==ParticleID);
    % control
    ParticleID = particle_vec_full(fall_ctrl_ids(i));
    pt_vec = [NaN(1,max(0,n_lags-fall_ctrl_ids(i)+1)) particle_vec_full(max(1,fall_ctrl_ids(i)-n_lags):fall_ctrl_ids(i))];
    dt_vec = [NaN(1,max(0,n_lags-fall_ctrl_ids(i)+1)) delta_vec_full(max(1,fall_ctrl_ids(i)-n_lags):fall_ctrl_ids(i))];
    fall_ctrl_array(i,pt_vec==ParticleID) = dt_vec(pt_vec==ParticleID);
end