% hmmm01_interpolate_data(project)
%
% DESCRIPTION
% Script to generate samples of local protein contentration at active gene
% loci and selected control locations
%
% ARGUMENTS
% project: master ID variable 
%
% OUTPUT: nucleus_struct: 

function nucleus_struct = hmmm01_interpolate_data(project,varargin)

dataPath = ['../dat/' project '/'];
minDp = 10;

for i = 1:numel(varargin)
    if strcmpi(varargin{i},'dropboxFolder')
        dataPath = [varargin{i+1} '\ProcessedEnrichmentData\' project '/'];
    end
    if ischar(varargin{i})
        if ismember(varargin{i},{'minDp','TresInterp'})
            eval([varargin{i} '=varargin{i+1};']);
        end
    end
end
load([dataPath '/nucleus_struct.mat'],'nucleus_struct')
% estimate interpolation time
med_time = nanmedian(diff([nucleus_struct.time]));
if exist('TresInterp')~= 1
    TresInterp = round(med_time);
end
interpGrid = 0:TresInterp:50*60;
%%% Cleaning Parameters
big_jump1 = prctile([nucleus_struct.fluo],99);
jump_threshold1 = big_jump1/1.5; % this should be a conservative threshold for single time step increase

for i = 1:length(nucleus_struct) 
    temp = nucleus_struct(i);
    trace1 = temp.fluo; %Load full trace, including intervening NaN's    
    pt_time = temp.time;      
    quality_flag = 1; % indicates whether trace suitable for inference
    
    if sum(~isnan(trace1)) == 0 
        t_start = 0;
        t_stop = -1;
    else
        t_start = interpGrid(find(interpGrid>=min(pt_time(~isnan(trace1))),1));
        t_stop = interpGrid(find(interpGrid<=max(pt_time(~isnan(trace1))),1,'last'));
    end
    time_interp = t_start:TresInterp:t_stop;
    
    if sum(~isnan(trace1)) < minDp
        quality_flag = 0;
        trace1_interp = NaN(size(time_interp));
        time_interp = NaN(size(time_interp));
    else
        %Null assumption is that all clusters of 6 or more NaNs are 0s. Smaller
        %clusters are assumed to have been missed nonzero dps
        trace1_nans = isnan(trace1);      
        %Look for clusters of 6 or more NaNs
        kernel = [1,1,1,1,1];
        tn_conv = conv(kernel,trace1_nans);
        tn_conv = tn_conv(3:end-2);
        z_ids = find(tn_conv==5);
        z_ids = unique([z_ids-1 z_ids z_ids+1]); % get set of z_ids    
        trace1(z_ids) = 0; % set clusters to zeros    
        trace1(trace1<0) = 0; % deal with negative values    
        % find single dp "blips". These will be replaced via interpolation
        tr_dd1 = abs([0 diff(diff(trace1)) 0]);
        trace1(tr_dd1>2*jump_threshold1) = NaN;    

        % interpolate remaining NaNs    
        query_points1 = pt_time(isnan(trace1));
        interp_t1 = pt_time(~isnan(trace1));
        interp_f1 = trace1(~isnan(trace1));
 
        new_f1 = interp1(interp_t1,interp_f1,query_points1);  
    
        trace1(ismember(pt_time,query_points1)) = new_f1;        

        %%% flag traces with unreasonably large rises or falls    
        tr_d1 = diff(trace1);    
        if max(abs(tr_d1)) >= jump_threshold1         
            quality_flag = 0;        
        end

        % Interpolate to standardize spacing    
        trace1_interp = interp1(pt_time,trace1,time_interp);    
    end
    interp_fields = {'xPos','yPos','ap_vector'};
    % interpolate other vector fields
    for j = 1:length(interp_fields)
        field_string = interp_fields{j};
        if isfield(nucleus_struct,field_string)
            init_vec = temp.(field_string);
            init_time = temp.time;   
            if numel(init_time) < 2 || numel(time_interp) < 2
                nucleus_struct(i).([interp_fields{j} '_interp']) = init_vec;                
            else
                nucleus_struct(i).([interp_fields{j} '_interp']) = interp1(init_time,init_vec,time_interp);
            end
        end
    end     
    if sum(isnan(trace1_interp)) > 0 && quality_flag == 1
        error('asfa')
    end
    nucleus_struct(i).fluo_interp = trace1_interp;    
    nucleus_struct(i).time_interp = time_interp;
    nucleus_struct(i).inference_flag = quality_flag;  
    nucleus_struct(i).TresInterp = TresInterp;
end
% save
save([dataPath '/nucleus_struct.mat'],'nucleus_struct')