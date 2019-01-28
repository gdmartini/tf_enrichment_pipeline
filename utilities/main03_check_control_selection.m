% main04_check_control_selection(project)
%
% DESCRIPTION
% Generates figures overlaying segmentation results onto histone channel
% and interface in which user can accept or reject spot assignments
%
% ARGUMENTS
% project: master ID variable

function main03_check_control_selection(project)

% specify paths
DataPath = ['../../dat/' project '/'];
SnipPath = [DataPath 'qc_images/'];
% load data
load(['../../dat/' project '/nucleus_struct_protein.mat']);
% snip_files = dir([SnipPath '*.mat']);
% check to see if nucleus structure already contains qc review info
if isfield(nucleus_struct_protein, 'qc_review_vec')
    qc_review_vec = [nucleus_struct_protein.qc_review_vec];
else 
    qc_review_vec = NaN(size([nucleus_struct_protein.xPos]));
end
qc_flag_vec = [nucleus_struct_protein.qc_flag_vec];
% set start frame
all_frames = find(~isnan(qc_flag_vec));
outstanding_frames = find(isnan(qc_review_vec)&(qc_flag_vec~=0&~isnan(qc_flag_vec)));
% generate indexing vectors
frame_index = [nucleus_struct_protein.frames];
set_index = [];
particle_index = [];
for i = 1:numel(nucleus_struct_protein)
    set_index = [set_index repelem(nucleus_struct_protein(i).setID, numel(nucleus_struct_protein(i).frames))];
    particle_index = [particle_index repelem(nucleus_struct_protein(i).ParticleID, numel(nucleus_struct_protein(i).frames))];
end    

% iterate through snip files
exit_flag = 0;
cm = jet(64);
index = outstanding_frames(1);
while ~exit_flag
    % create sister_struct(i) struct    
    frame = frame_index(index);
    ParticleID = particle_index(index);
    setID = set_index(index);    
    % load snip data
    load([SnipPath 'pt' num2str(1e4*ParticleID) '_frame' sprintf('%03d',frame) '.mat']);
    
    %%%%%%%%%%%%%%%%%%%%%%% load image stack %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    cc = '';
    while ~strcmp(cc,'0')&&~strcmp(cc,'1')      
        dist_snip = qc_spot.dist_snip;
        dist_rescaled = 1 + ceil(dist_snip/max(dist_snip(:)) * 63);
        cm = jet(64);
        cm = [[1 1 1] ; cm(2:end,:)];
        dist_rgb = ind2rgb(dist_rescaled,cm);
        % get frame center
        x_center = qc_spot.x_center;
        y_center = qc_spot.y_center;
        
        qc_fig = figure;                 
        imshow(imadjust(mat2gray(qc_spot.mcp_snip)),'InitialMagnification','fit');                        
        hold on
        p = imshow(dist_rgb);        
        p.AlphaData = .4;        
        scatter(qc_spot.xp-x_center,qc_spot.yp-y_center,30,'MarkerFaceColor',cm(30,:),'MarkerEdgeAlpha',0)
        scatter(qc_spot.xc-x_center,qc_spot.yc-y_center,30,'MarkerFaceColor',cm(60,:),'MarkerEdgeAlpha',0)
          
        ct=waitforbuttonpress;
        cc=get(SpotFig,'currentcharacter');
        if strcmp(cc,'1')||strcmp(cc,'0')           
            particle_labels(i) = eval(cc);
            iter = min(iter+1,numel(indices));           
        elseif strcmp(cc,'x')
            exit_flag = 1;
            break
        elseif strcmp(cc,'a')
            nc = min(zDim,nc+1);
        elseif strcmp(cc,'z')
            nc = max(1,nc-1);      
        elseif strcmp(cc,'n')
            iter = max(1,iter-1);
            break
        elseif strcmp(cc,'m')
            iter = min(numel(indices),iter+1);
            break
        end       
    end 
    close all
    if exit_flag
        disp('Exiting')
        break
    end
    i = indices(iter);
end
