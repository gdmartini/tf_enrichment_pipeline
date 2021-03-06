% Script to validate methods for discriminating between target and
% off-target loci using protein concentration
clear
close all
% set id variables and load data
dropboxFolder = 'E:\Nick\Dropbox (Garcia Lab)\ProcessedEnrichmentData\';
project_cell = {'Dl_Venus_snaBAC_mCherry_LeicaZoom2x_7uW14uW_02','Dl_Venus_hbP2P_mCherry_LeicaZoom2x_HighPower'};
master_struct = struct;
for i = 1:numel(project_cell)
    master_struct(i).project = project_cell{i};
    load([dropboxFolder project_cell{i} '/nucleus_struct_protein.mat'])
    master_struct(i).nucleus_struct_protein = nucleus_struct_protein;
    
    % get other project info
    project = project_cell{i};
    underscores = strfind(project,'_');
    master_struct(i).protein_name = project(1:underscores(1)-1);
    master_struct(i).protein_fluor = project(underscores(1)+1:underscores(2)-1);
    master_struct(i).gene_name = project(underscores(2)+1:underscores(3)-1);
    if numel(underscores) == 3
        ind = numel(project);
    else
        ind = underscores(4)-1;
    end
    master_struct(i).gene_fluor = project(underscores(3)+1:ind);
end


% assign pseudo sister pairings from across two sets
mindp = 20; % only nuclei with 20 time-step temporal overlap will be permitted
% ugh need to interpolate...
interpGrid = 0:20:50*60;
interp_fields = {'fluo','mf_null_protein_vec','spot_protein_vec'};
for i = 1:numel(master_struct)
    nucleus_struct_protein = master_struct(i).nucleus_struct_protein;
    for j = 1:numel(nucleus_struct_protein)
        nc_time = nucleus_struct_protein(j).time;
        time_interp = interpGrid(interpGrid<=nc_time(end)&interpGrid>=nc_time(1));
        nucleus_struct_protein(j).time_interp = time_interp;
        qc_flag = true;
        if numel(nc_time) > 1
            for k = 1:numel(interp_fields)
                vec = nucleus_struct_protein(j).(interp_fields{k});
                vec_ft = ~isnan(vec);
                qc_flag = qc_flag & sum(vec_ft) > mindp;
                if sum(vec_ft) > 1 
                    nucleus_struct_protein(j).([interp_fields{k} '_interp']) = ...
                        interp1(nc_time(vec_ft),vec(vec_ft),time_interp);

                else
                    nucleus_struct_protein(j).([interp_fields{k} '_interp']) = NaN(size(time_interp));
                end
            end
            nucleus_struct_protein(j).time_interp = time_interp;
        else
            qc_flag = false;
            for k = 1:numel(interp_fields)
                vec = nucleus_struct_protein(j).(interp_fields{k});                
                nucleus_struct_protein(j).([interp_fields{k} '_interp']) = vec;
            end
            time_interp = nc_time;
            nucleus_struct_protein(j).time_interp = time_interp;
        end
        nucleus_struct_protein(j).qc_flag = qc_flag;
        nucleus_struct_protein(j).first_time = time_interp(1);
        nucleus_struct_protein(j).last_time = time_interp(end);
    end
    master_struct(i).nucleus_struct_protein = nucleus_struct_protein;    
end
% generate index vectors of nuclei in each set that are up to snuff
id_cell = cell(size(project_cell));
dims = NaN(size(project_cell));
for i = 1:numel(project_cell)
    qc_vec = [master_struct(i).nucleus_struct_protein.qc_flag];
    id_cell{i} = find(qc_vec);
    dims(i) = numel(id_cell{i});
end
start_array = NaN(dims(1),dims(2),2);
stop_array = NaN(dims(1),dims(2),2);
options = [1 2];
for i = 1:numel(master_struct)
    start_vec = [master_struct(i).nucleus_struct_protein(id_cell{i}).first_time];    
    stop_vec = [master_struct(i).nucleus_struct_protein(id_cell{i}).last_time];
    if i == 1
        start_array(:,:,i) = repmat(start_vec,dims(options(options~=i)),1)'; % wil probs throw error
        stop_array(:,:,i) = repmat(stop_vec,dims(options(options~=i)),1)'; % wil probs throw error
    else
        start_array(:,:,i) = repmat(start_vec,dims(options(options~=i)),1); % wil probs throw error
        stop_array(:,:,i) = repmat(stop_vec,dims(options(options~=i)),1); % wil probs throw error
    end
end
% calculate overlaps
overlap_mat = min(stop_array,[],3) - max(start_array,[],3);
overlap_flag_array = overlap_mat >= mindp*20;
% calculate average per-point difference between mf protein vectors for
% those nuclei that overlap sufficently often

[id_vec1,id_vec2] = find(overlap_flag_array);
delta_array = zeros(size(overlap_mat));
for i = 1:numel(id_vec1) 
    t1 = master_struct(1).nucleus_struct_protein(id_cell{1}(id_vec1(i))).time_interp;
    t2 = master_struct(2).nucleus_struct_protein(id_cell{2}(id_vec2(i))).time_interp;
    t1_filter = ismember(t1,t2);
    t2_filter = ismember(t2,t1);
    mf_delta_vec = master_struct(1).nucleus_struct_protein(id_cell{1}(id_vec1(i))).mf_null_protein_vec_interp(t1_filter) - ...
        master_struct(2).nucleus_struct_protein(id_cell{2}(id_vec2(i))).mf_null_protein_vec_interp(t2_filter);    
    if sum(~isnan(mf_delta_vec)) > mindp
        delta_array(id_vec1(i),id_vec2(i)) = nanmean(abs(mf_delta_vec));
    end     
end

% use existing code to solve linear assignment problem
M = matchpairs(delta_array, 1);
% SICK. Now, finally, calculate error rate for pairs using simple binary
% decision procedure (higher=target)
dev_vec = NaN(1,size(M,1));
success_mat = NaN(size(M,1),round(nanmax(overlap_mat(:))/20));
for i = 1:size(M,1)
    t1 = master_struct(1).nucleus_struct_protein(id_cell{1}(M(i,1))).time_interp;
    t2 = master_struct(2).nucleus_struct_protein(id_cell{2}(M(i,2))).time_interp;
    t1_filter = ismember(t1,t2);
    t2_filter = ismember(t2,t1);
    % extract mf protein vectors
    mf1 = master_struct(1).nucleus_struct_protein(id_cell{1}(M(i,1))).mf_null_protein_vec_interp(t1_filter);
    mf2 = master_struct(2).nucleus_struct_protein(id_cell{2}(M(i,2))).mf_null_protein_vec_interp(t2_filter);
    % extract spot protein vectors
    sp1 = master_struct(1).nucleus_struct_protein(id_cell{1}(M(i,1))).spot_protein_vec_interp(t1_filter);
    sp2 = master_struct(2).nucleus_struct_protein(id_cell{2}(M(i,2))).spot_protein_vec_interp(t2_filter);
    % take intersection
    pt_ft = ~isnan(mf1) & ~isnan(mf2) & ~isnan(sp1) & ~isnan(sp2);    
    target_delta = sp1(pt_ft) - mf1(pt_ft);
    control_delta = sp2(pt_ft) - mf2(pt_ft);
    n_vec = 1:numel(control_delta);
    diff_vec = (cumsum(target_delta) - cumsum(control_delta)) ./ n_vec;
    success_vec = diff_vec > 0;
    success_mat(i,1:numel(success_vec)) = success_vec;
end
%%
