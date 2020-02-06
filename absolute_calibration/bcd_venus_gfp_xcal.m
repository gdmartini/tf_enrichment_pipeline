clear
close all

addpath('../utilities/')
% define basic path variables
DropboxFolder =  'E:\Nick\LivemRNA\Dropbox (Personal)\';
Bcd_GFP_project = 'Bcd-GFP_hbP2P-mCh';
Bcd_Venus_project = 'Bcd-Venus';


[RawResultsRoot, ~, ~] =   header_function(DropboxFolder, Bcd_GFP_project);
DataPath =  [DropboxFolder 'ProcessedEnrichmentData\absolute_calibration\'];
mkdir(DataPath);
FigPath = [DropboxFolder 'LocalEnrichmentFigures/PipelineOutput/absolute_calibration/venus_gfp_xcal/'];
mkdir(FigPath);

%% load absolute Bcd concentration data from Gregor 2007
bcd_abs_path = 'E:\Nick\LivemRNA\Dropbox (Personal)\processedenrichmentdata\absolute_calibration\GregorData2007\';
bkg_data = readtable([bcd_abs_path 'GregorBkg2007.csv']);
bkg_data_clean.AP = bkg_data.Var1;
bkg_data_clean.nM = bkg_data.Var2;
bcd_data = readtable([bcd_abs_path 'GregorBcd2007.csv']);
bcd_data_clean.AP = bcd_data.Var1(1:41);
bcd_data_clean.nM = bcd_data.Var2(1:41);
%%
% find sheets
sheet_path = [RawResultsRoot 'DataStatus.xlsx'];
[~,sheet_names]=xlsfinfo(sheet_path);

sheet_index_gfp = find(ismember(sheet_names,Bcd_GFP_project));
sheet_index_venus = find(ismember(sheet_names,Bcd_Venus_project));

% get prefix names
% GFP
[~,~,sheet_cell] = xlsread(sheet_path,sheet_index_gfp);
name_col = sheet_cell(1:33,1); % hard coded for now
ready_ft = contains(name_col,'ReadyForEnrichment');
ready_cols = 1 + find([sheet_cell{ready_ft,2:end}]==1);
sheet_cell = sheet_cell(:,[1 ready_cols]);
% get list of project names
prefix_ft = contains(name_col,'Prefix');
prefix_cell_raw = sheet_cell(prefix_ft,2:end);
prefix_cell_gfp = {};
for i = 1:numel(prefix_cell_raw)
    if ~isempty(prefix_cell_raw{i})
        eval([prefix_cell_raw{i} ';'])
        prefix_cell_gfp = [prefix_cell_gfp{:} {Prefix}];
    end
end

% Venus
[~,~,sheet_cell] = xlsread(sheet_path,sheet_index_venus);
name_col = sheet_cell(1:33,1); % hard coded for now
ready_ft = contains(name_col,'ReadyForEnrichment');
ready_cols = 1 + find([sheet_cell{ready_ft,2:end}]==1);
sheet_cell = sheet_cell(:,[1 ready_cols]);
% get list of project names
prefix_ft = contains(name_col,'Prefix');
prefix_cell_raw = sheet_cell(prefix_ft,2:end);
prefix_cell_venus = {};
for i = 1:numel(prefix_cell_raw)
    if ~isempty(prefix_cell_raw{i})
        eval([prefix_cell_raw{i} ';'])
        prefix_cell_venus = [prefix_cell_venus{:} {Prefix}];
    end
end
% get pixel size (should be consistent across sets
load([RawResultsRoot prefix_cell_venus{1} '/FrameInfo.mat']) 
PixelSize = FrameInfo(1).PixelSize;
IntegrationRadius=2;       %Radius of the integration region in um
IntegrationRadius=floor(IntegrationRadius/FrameInfo(1).PixelSize); %Radius of the integration in pixels
if ~mod(IntegrationRadius,2)
    IntegrationRadius=IntegrationRadius+1;
end
Circle=false(3*IntegrationRadius,3*IntegrationRadius);
Circle=MidpointCircle(Circle,IntegrationRadius,1.5*IntegrationRadius+0.5,...
    1.5*IntegrationRadius+0.5,1);
n_pixels = sum(Circle(:));

% load nuclear fluorescence data

% Venus
time_vec_venus = [];
ap_vec_venus = [];
pt_vec_venus = [];
set_vec_venus = [];

for i = 1:numel(prefix_cell_venus)
    % load
    load([RawResultsRoot prefix_cell_venus{i} '/CompiledNuclei.mat']) 
    % extract data
    nc14_time = ElapsedTime(nc14);
    TimeMat = repmat(ElapsedTime',1,size(AllTracesVector,2))-nc14_time;
    APMat = repmat(AllTracesAP',size(AllTracesVector,1),1);
    % generate vectors
    qc_filter = ~isnan(AllTracesVector)&TimeMat>=0;
    time_vec = TimeMat(qc_filter)';
    ap_vec = APMat(qc_filter)';
    pt_vec = AllTracesVector(qc_filter)';
    % add to master vectors
    time_vec_venus = [time_vec_venus time_vec];
    ap_vec_venus = [ap_vec_venus ap_vec];
    pt_vec_venus = [pt_vec_venus pt_vec/n_pixels];
    set_vec_venus = [set_vec_venus repelem(i,numel(time_vec))];
end

% Bcd
time_vec_gfp = [];
ap_vec_gfp = [];
pt_vec_gfp = [];
set_vec_gfp = [];

for i = 1:numel(prefix_cell_gfp)
    % load
    load([RawResultsRoot prefix_cell_gfp{i} '/CompiledNuclei.mat']) 
    % extract data
    nc14_time = ElapsedTime(nc14);
    TimeMat = repmat(ElapsedTime',1,size(AllTracesVector,2))-nc14_time;
    APMat = repmat(AllTracesAP',size(AllTracesVector,1),1);
    % generate vectors
    qc_filter = ~isnan(AllTracesVector)&TimeMat>=0;
    time_vec = TimeMat(qc_filter)';
    ap_vec = APMat(qc_filter)';
    pt_vec = AllTracesVector(qc_filter)';
    % add to master vectors
    time_vec_gfp = [time_vec_gfp time_vec];
    ap_vec_gfp = [ap_vec_gfp ap_vec];
    pt_vec_gfp = [pt_vec_gfp pt_vec/n_pixels];
    set_vec_gfp = [set_vec_gfp repelem(i,numel(time_vec))];
end

% Make basic plots for each fluor
symbol_cell = {'o','s','x','*'};
close all

% Venus 
bcd_venus_fig = figure;
cmap = flipud(brewermap([],'Spectral'));
colormap(cmap);
hold on
% for i = 1:numel(prefix_cell_venus)
%     sft = set_vec_venus == i;
    scatter(ap_vec_venus,pt_vec_venus,15,time_vec_venus,'o','filled','MarkerEdgeAlpha',0,'MarkerFaceAlpha',.2);
% end
h = colorbar;
xlabel('% embryo length')
ylabel('Bcd-Venus intensity (au)')
ylabel(h,'minutes into nc14')
xlim([.17 .38])
set(gca,'FontSize',14)
saveas(bcd_venus_fig,[FigPath 'bcd_venus_scatter.png'])

% GFP
bcd_gfp_fig = figure;
colormap(cmap);
hold on
% for i = 1:numel(prefix_cell_gfp)
%     sft = set_vec_gfp == i;
    scatter(ap_vec_gfp,pt_vec_gfp,15,time_vec_gfp,'filled','MarkerEdgeAlpha',0,'MarkerFaceAlpha',.2);
% end
h = colorbar;
xlabel('% embryo length')
ylabel('Bcd-GFP intensity (au)')
ylabel(h,'minutes into nc14')
% xlim([.17 .38])
set(gca,'FontSize',14)
saveas(bcd_gfp_fig,[FigPath 'bcd_gfp_scatter.png'])

% Perform simple cross-calibration using average values
ap_res = 1;
time_res = 1;
ap_range = 15:50;
time_range = 5:25;
time_mat = repmat(time_range',1,numel(ap_range));

venus_fluo_mat = NaN(numel(time_range),numel(ap_range));
venus_ct_mat = NaN(numel(time_range),numel(ap_range));
gfp_fluo_mat = NaN(numel(time_range),numel(ap_range));
gfp_ct_mat = NaN(numel(time_range),numel(ap_range));

for  a = 1:numel(ap_range)
    for t = 1:numel(time_range)
        v_filter = round(100*ap_vec_venus)==ap_range(a) & round(time_vec_venus)==time_range(t);
        venus_fluo_mat(t,a) = nanmean(pt_vec_venus(v_filter));
        venus_ct_mat(t,a) = sum(v_filter);
        
        g_filter = round(100*ap_vec_gfp)==ap_range(a) & round(time_vec_gfp)==time_range(t);
        gfp_fluo_mat(t,a) = nanmean(pt_vec_gfp(g_filter));
        gfp_ct_mat(t,a) = sum(g_filter);
    end
end
%%
close all
xcal_fig = figure;
colormap(cmap)
scatter(venus_fluo_mat(:),gfp_fluo_mat(:),20,time_mat(:),'filled')
grid on 
h = colorbar;
xlabel('Bcd-Venus intensity (au)')
ylabel('Bcd-GFP intensity (au)')
ylabel(h,'minutes into nc14')
set(gca,'FontSize',14)
axis([0 .4 0 .4])
saveas(xcal_fig,[FigPath 'venus_gfp_xcal_scatter.png'])


% 

%% Compare to Gregor Data

% obtain estimate as function of AP
ap_vec_bcd = 100*bcd_data_clean.AP;
ap_vec_bkg = 100*bkg_data_clean.AP;
bcd_nM_raw = interp1(ap_vec_bcd,bcd_data_clean.nM,ap_range);
bkg_nM_raw = interp1(ap_vec_bkg,bkg_data_clean.nM,ap_range);

bcd_nM_vec = bcd_nM_raw - bkg_nM_raw;
% pull Bcd-Venus data at t=15
bcd_venus_vec = venus_fluo_mat(time_range==15,:);

md = fitlm(bcd_nM_vec, bcd_venus_vec);

%% calculate confocal volume 
rxy = .61*510e-9/1.4;
rz = 2*1.5*510e-9/1.4;
% kappa = 2.33 * 1.5 / 1.4;
VoxelSize = PixelSize^2 * .5;
micronLiter = 1e-15;
Nmol = 6.022e23;

20 * Nmol * 1e-9 * VoxelSize / 1e15 

