% main02_sample_local_protein(projectName, varargin)
%
% DESCRIPTION
% Script to generate samples of local protein contentration at active gene
% loci and selected control locations
%
% ARGUMENTS
% projectName: master ID variable (should match a tab name in the 
%              DataStatus.xlsx spreadsheet)
%
% OPTIONS
% script allows any default variable to be set using format:
%       "VariableNameString", VariableValue
%
% OUTPUT: spot_struct_protein: compiled data set with protein samples

function spot_struct_protein = main02_sample_local_protein(projectName,varargin)

    %% %%%%%%%%%%%%%%%%%%%%%%% Set Defaults %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    addpath(genpath('utilities'))
    close all force

    ROIRadiusSpot_um = .2; % radus (um) of region used to query and compare TF concentrations
    minSampleSep_um = 1.5; %um
    minEdgeSep_um = .25; %um
    snippet_size_um = 1.5;
    min_nucleus_radius_um = 2;
    max_nucleus_radius_um = 4;

    segmentNuclei = 0;
    use3DSpotInfo = 1;
    NumWorkers = 24;
    % PSF info for 3D sampling
    use_psf_fit_dims = false; % NL: currently no supported
    xy_sigma_um = 0.25;% um 
    xy_sigma_nuclear_um = 1.5;
    z_sigma_um = 0.6; % um
    ignoreQC = false;
    write_snip_flag = false; %NL: what does this do?

    %% %%%%%%%%%%%%%%%%%%%%%%% Check for optional inputs %%%%%%%%%%%%%%%%%%%%%%
    for i = 1:(numel(varargin)-1)  
        if i ~= numel(varargin)        
            eval([varargin{i} '=varargin{i+1};']);                
        end    
    end

    %% %%%%%%%%%%%%%%%%%%%%%%% Save key sampling parameters %%%%%%%%%%%%%%%%%%%
    proteinSamplingInfo = struct;
    proteinSamplingInfo.ROIRadiusSpot = ROIRadiusSpot_um;
    proteinSamplingInfo.minSampleSep_um = minSampleSep_um;
    proteinSamplingInfo.minEdgeSep_um = minEdgeSep_um; 
    proteinSamplingInfo.xy_sigma_nuclear_um = xy_sigma_nuclear_um;
    proteinSamplingInfo.ROIRadiusSpot_um = ROIRadiusSpot_um;
    proteinSamplingInfo.snippet_size_um = snippet_size_um;
    proteinSamplingInfo.xy_sigma_um = xy_sigma_um;
    proteinSamplingInfo.z_sigma_um = z_sigma_um;
    proteinSamplingInfo.min_nucleus_radius_um = min_nucleus_radius_um;
    proteinSamplingInfo.max_nucleus_radius_um = max_nucleus_radius_um;

    %% %%%%%%%%%%%%%%%%%%%%%%% Get project info %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [liveProject, ~, nucleusName, hasAPInfo, has3DSpotInfo] = headerFunction(projectName);
    use3DSpotInfo = use3DSpotInfo&&has3DSpotInfo;
    proteinSamplingInfo.use3DSpotInfo = use3DSpotInfo;

    %% %%%%%%%%%%%%%%%%%%%%%%% Load data and clean trace data %%%%%%%%%%%%%%%%%
    load(nucleusName,'nucleus_struct')

    if false % NL: Currently not working
        load([liveProject.dataPath '/psf_dims.mat'],'psf_dims')
    end

    % make paths 
    snipPath = [liveProject.dataPath '/qc_images/'];
    refPath = [liveProject.dataPath '/refFrames/'];
    mkdir(refPath)
    mkdir(snipPath)

    % remove frames where no particle was observed
    spot_struct_protein = truncateParticleFields(nucleus_struct,use3DSpotInfo,hasAPInfo);

    %% %%%%%%%%%%%%%%%%%%%%%%% Generate indexing vectors  %%%%%%%%%%%%%%%%%%%%%
    [RefStruct, set_frame_array] = generateReferenceVectors(spot_struct_protein,refPath,use3DSpotInfo,ignoreQC);

    qc_structure = struct;


    %% %%%%%%%%%%%%%%%%%%%%%%% Initialize enrichment-related fields  %%%%%%%%%%

    [spot_struct_protein, ~] = initializeProteinFields(spot_struct_protein, use3DSpotInfo);

    snip_fields = {'spot_protein_snips', 'edge_control_protein_snips',...
        'spot_mcp_snips','edge_control_mcp_snips'};


    %% %%%%%%%%%%%%%%%%%%%%%%% Nucleus segmentation  %%%%%%%%%%%%%%%%%%%%%%%%%%

    % first check to see if segmentation files exist
    [segmentNuclei, segmentIndices] = ...
              handleSegmentationOptions(RefStruct,segmentNuclei);

    if segmentNuclei
        disp('segmenting nuclei...')    
        nuclearSegmentation(liveProject, RefStruct, segmentIndices, NumWorkers);      
    end

    %% %%%%%%%%%%%%%%%%%%%%%%% Local Protein Sampling %%%%%%%%%%%%%%%%%%%%%%%%%
    
    % try turning the structure into a cell array to make it more ammenable
    % to parallelization
    spot_cell = cell(1,length(spot_struct_protein));
    for i = 1:length(spot_cell)
        spot_temp = spot_struct_protein(i);
        spot_cell{i} = spot_temp;
    end
        
    % convert structure to cell array for conveneince
    pool = gcp('nocreate');
    if isempty(pool)
      parpool(12)
    end

    D = parallel.pool.DataQueue;    
    afterEach(D, @nUpdateWaitbar);

%     N = size(set_frame_array,1);
%     p = 1;
% 
%     h = waitbar(0,'Sampling local protein...');
    for i = 1:size(set_frame_array,1)            
    %     waitbar(i/size(set_frame_array,1),h)

        samplingInfo = struct; % use this structure to keep track of key tracking-related info

        % read basic info from set_frame array
        samplingInfo.SetID = set_frame_array(i,1);
        samplingInfo.Frame = set_frame_array(i,2);  

        % get channel info
        Prefix = liveProject.includedExperimentNames{samplingInfo.SetID}; 
        currExperiment = liveProject.includedExperiments{samplingInfo.SetID};                  
        samplingInfo.proteinChannel = currExperiment.inputChannels;
        samplingInfo.mcpChannel = currExperiment.spotChannels;

        %% %%%%%%%%%%%%%%%%%%%%% Set size parameters  %%%%%%%%%%%%%%%%%%%%%%%%%
        currExperiment = LiveExperiment(Prefix);
        PixelSize = currExperiment.pixelSize_nm / 1e3;  
        zStep = currExperiment.zStep_um;    

        if length(samplingInfo.mcpChannel) > 1 
          error('This pipeline does not currently support multiple spot channels')
        elseif length(samplingInfo.proteinChannel) > 1
          error('This pipeline does not currently support multiple protein channels')
        end

        % Generate reference vectors
        xDim = currExperiment.xDim;
        yDim = currExperiment.yDim;
        zDim = currExperiment.zDim;
        [samplingInfo.x_ref,samplingInfo.y_ref,samplingInfo.z_ref] = meshgrid(1:xDim,1:yDim,1:zDim);

        % calculate basic parameters for sampling
        sampParamNames = fieldnames(proteinSamplingInfo);   
        for s = 1:length(sampParamNames)
          paramName = sampParamNames{s};
          samplingInfo.(paramName(1:end-3)) = proteinSamplingInfo.(paramName) / PixelSize;
        end
        samplingInfo.z_sigma = proteinSamplingInfo.z_sigma_um / zStep;
        samplingInfo.min_nucleus_area = pi*samplingInfo.min_nucleus_radius^2;
        samplingInfo.max_nucleus_area = pi*samplingInfo.max_nucleus_radius^2;
        samplingInfo.snippet_size = round(samplingInfo.snippet_size);

        % calculate characteristic drift to use for simulated spot
        samplingInfo.driftTol = calculateVirtualSpotDrift(RefStruct,PixelSize);

        % load spot and nucleus reference frames
        nc_ref_name = [refPath 'nc_ref_frame_set' sprintf('%02d',samplingInfo.SetID) '_frame' sprintf('%03d',samplingInfo.Frame) '.mat'];
        temp = load(nc_ref_name,'nc_ref_frame');
        nc_ref_frame = temp.nc_ref_frame;        
        samplingInfo.nc_dist_frame = bwdist(~nc_ref_frame);    
        
        spot_ref_name = [refPath 'spot_roi_frame_set' sprintf('%02d',samplingInfo.SetID) '_frame' sprintf('%03d',samplingInfo.Frame) '.mat'];
        temp = load(spot_ref_name,'spot_dist_frame');    
        spot_dist_frame = temp.spot_dist_frame;
        samplingInfo.spot_dist_frame = spot_dist_frame;

        % get indices of particles in current set/frame 
        frame_set_filter = RefStruct.setID_ref==samplingInfo.SetID&RefStruct.frame_ref==samplingInfo.Frame;
        samplingInfo.frame_set_indices = find(frame_set_filter);                

        % load stacks    
        proteinPath = [currExperiment.preFolder  Prefix '_' sprintf('%03d',samplingInfo.Frame) '_ch0' num2str(samplingInfo.proteinChannel) '.tif'];
        samplingInfo.protein_stack = imreadStack(proteinPath);
        if size(samplingInfo.protein_stack,3) == zDim+2
          samplingInfo.protein_stack = samplingInfo.protein_stack(:,:,2:end); 
        else
          error('Unrecognized form of z-padding')
        end    


        mcpPath = [currExperiment.preFolder  Prefix '_' sprintf('%03d',samplingInfo.Frame) '_ch0' num2str(samplingInfo.mcpChannel) '.tif'];
        samplingInfo.mcp_stack = imreadStack2(mcpPath);
        if size(samplingInfo.mcp_stack,3) == zDim+2
          samplingInfo.mcp_stack = samplingInfo.mcp_stack(:,:,2:end);
        else
          error('Unrecognized form of z-padding')
        end

        % generate lookup table of inter-nucleus distances
        nc_x_vec = RefStruct.nc_x_ref(frame_set_filter);
        nc_y_vec = RefStruct.nc_y_ref(frame_set_filter);  
        x_dist_mat = repmat(nc_x_vec,numel(nc_x_vec),1)-repmat(nc_x_vec',1,numel(nc_x_vec));
        y_dist_mat = repmat(nc_y_vec,numel(nc_y_vec),1)-repmat(nc_y_vec',1,numel(nc_y_vec));
        samplingInfo.r_dist_mat = sqrt(double(x_dist_mat).^2 + double(y_dist_mat).^2);            

        % initialize temporary arrays to store snip info 
        temp_snip_struct = struct;
        for j = 1:length(snip_fields)
            temp_snip_struct.(snip_fields{j})  = NaN(2*samplingInfo.snippet_size+1,2*samplingInfo.snippet_size+1,length(nc_x_vec));
        end    
        
        % iterate through spots
        qc_mat = struct;
        j_pass = 1;
        for j = samplingInfo.frame_set_indices        
            % get indexing info         
            samplingInfo.spotIndex = RefStruct.particle_index_ref(j);
            samplingInfo.spotSubIndex = RefStruct.particle_subindex_ref(j);

            % get location info
            x_nucleus = RefStruct.nc_x_ref(j);
            y_nucleus = RefStruct.nc_y_ref(j);

            x_spot = RefStruct.spot_x_ref(j);
            x_index = min([ xDim max([1 round(x_spot)])]);
            y_spot = RefStruct.spot_y_ref(j);
            y_index = min([ yDim max([1 round(y_spot)])]);
            z_spot = RefStruct.spot_z_ref(j)-1.0; % adjust for z padding

            % extract mask 
            nucleus_mask = nc_ref_frame == RefStruct.master_nucleusID_ref(j);             

            %% %%%%%%%%%%%%%%%%%%% Sample protein levels %%%%%%%%%%%%%%%%%%%%%%

            % sample protein near locus
            spot_cell{samplingInfo.spotIndex}.spot_protein_vec(samplingInfo.spotSubIndex) = ...
                  sample_protein_3D(samplingInfo,samplingInfo.protein_stack,...
                  x_spot,y_spot,z_spot,samplingInfo.xy_sigma,samplingInfo.z_sigma);

            spot_struct_protein(samplingInfo.spotIndex).spot_mcp_vec(samplingInfo.spotSubIndex) = ...
                  sample_protein_3D(samplingInfo,samplingInfo.mcp_stack,...
                  x_spot,y_spot,z_spot,samplingInfo.xy_sigma,samplingInfo.z_sigma);

            % make sure size is reasonable and that spot is inside nucleus
            if sum(nucleus_mask(:)) < samplingInfo.min_nucleus_area || sum(nucleus_mask(:)) > samplingInfo.max_nucleus_area...
                || ~nucleus_mask(y_index,x_index)

                spot_struct_protein(samplingInfo.spotIndex).edge_qc_flag_vec(samplingInfo.spotSubIndex) = -1;            
                spot_struct_protein(samplingInfo.spotIndex).serial_qc_flag_vec(samplingInfo.spotSubIndex) = -1;

                continue
            end 

            % sample snippets        
            temp_snip_struct.spot_protein_snips(:,:,j_pass) = sample_snip_3D(x_spot,y_spot,z_spot,samplingInfo,samplingInfo.protein_stack);
            temp_snip_struct.spot_mcp_snips(:,:,j_pass) = sample_snip_3D(x_spot,y_spot,z_spot,samplingInfo,samplingInfo.mcp_stack); 

            % Take average across all pixels within 1.5um of nuclues center           
            spot_struct_protein(samplingInfo.spotIndex).nucleus_protein_vec(samplingInfo.spotSubIndex) = sample_protein_3D(...
              samplingInfo,samplingInfo.protein_stack,x_nucleus,y_nucleus,z_spot,samplingInfo.xy_sigma_nuclear,samplingInfo.z_sigma);        

            %% %%%%%%%%%%%% Draw edge control spot %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            spot_struct_protein = findEdgeControlWrapper(spot_struct_protein,...
              RefStruct,samplingInfo,nucleus_mask,x_index,y_index,samplingInfo.spotIndex,samplingInfo.spotSubIndex,j_pass);

            % Draw control samples (as appropriate)    
            if spot_struct_protein(samplingInfo.spotIndex).edge_qc_flag_vec(samplingInfo.spotSubIndex) > 0  
                edge_control_x = spot_struct_protein(samplingInfo.spotIndex).edge_null_x_vec(samplingInfo.spotSubIndex);
                edge_control_y = spot_struct_protein(samplingInfo.spotIndex).edge_null_y_vec(samplingInfo.spotSubIndex);     

                spot_struct_protein(samplingInfo.spotIndex).edge_null_protein_vec(samplingInfo.spotSubIndex) = sample_protein_3D(...
                  samplingInfo,samplingInfo.protein_stack,edge_control_x,edge_control_y,z_spot,...
                  samplingInfo.xy_sigma,samplingInfo.z_sigma);

                spot_struct_protein(samplingInfo.spotIndex).edge_null_mcp_vec(samplingInfo.spotSubIndex) = sample_protein_3D(...              
                  samplingInfo,samplingInfo.mcp_stack,edge_control_x,edge_control_y,z_spot,...
                  samplingInfo.xy_sigma,samplingInfo.z_sigma);

                % draw snips    
                temp_snip_struct.edge_control_protein_snips(:,:,j_pass) = sample_snip_3D(edge_control_x,edge_control_y,z_spot,samplingInfo,samplingInfo.protein_stack);
                temp_snip_struct.edge_control_mcp_snips(:,:,j_pass) = sample_snip_3D(edge_control_x,edge_control_y,z_spot,samplingInfo,samplingInfo.mcp_stack);
            end                  

            %% %%%%%%%%%%%% Draw serialized control spot %%%%%%%%%%%%%%%%%%%%%%%           
            spot_struct_protein = drawSerializedControlSpot(...
              spot_struct_protein,samplingInfo,x_index,y_index,samplingInfo.spotIndex,samplingInfo.spotSubIndex,nucleus_mask...
              );
            
            % sample protein 
            serial_control_x = spot_struct_protein(samplingInfo.spotIndex).serial_null_x_vec(samplingInfo.spotSubIndex);
            serial_control_y = spot_struct_protein(samplingInfo.spotIndex).serial_null_y_vec(samplingInfo.spotSubIndex);

            spot_struct_protein(samplingInfo.spotIndex).serial_null_protein_vec(samplingInfo.spotSubIndex) = sample_protein_3D(...
                  samplingInfo,samplingInfo.protein_stack,serial_control_x,serial_control_y,z_spot,samplingInfo.xy_sigma,samplingInfo.z_sigma);
            spot_struct_protein(samplingInfo.spotIndex).serial_null_mcp_vec(samplingInfo.spotSubIndex) = sample_protein_3D(...
                  samplingInfo,samplingInfo.mcp_stack,serial_control_x,serial_control_y,z_spot,samplingInfo.xy_sigma,samplingInfo.z_sigma);


            %% %%%%%%%%%%%% Check for sister spot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %         x_spot_sister = NaN;
    %         y_spot_sister = NaN;
    %         nucleus_indices = find(refVecStruct.master_nucleusID_ref(j)==refVecStruct.master_nucleusID_ref(frame_set_filter));
    %         if length(nucleus_indices) == 2            
    %             sister_index = nucleus_indices(nucleus_indices~=j_pass);
    %             x_spot_sister = spot_x_vec(sister_index);
    %             y_spot_sister = spot_y_vec(sister_index);
    %         end            

            %% %%%%%%%%%%%%% save qc data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                 
    %         qc_mat = updateQCmat;

            % increment
            j_pass = j_pass + 1;        
        end 
        qc_structure(i).qc_mat = fliplr(qc_mat);  

        %% %%%%%%%%%%%%%%%%%%%% save snip data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

        % initialize struct to store snip data
        snip_data = struct;        
        % add protein snip stacks    
        for j = 1:length(snip_fields)
            snip_data.(snip_fields{j}) = temp_snip_struct.(snip_fields{j});
        end    

        % store key ID variables
        snip_data.frame = samplingInfo.Frame;
        snip_data.setID = samplingInfo.SetID;
        snip_data.particle_id_vec = RefStruct.particleID_ref(frame_set_filter);
    
        % indexing vectors    
        snip_data.nc_sub_index_vec = RefStruct.particle_subindex_ref(frame_set_filter);
        snip_data.nc_lin_index_vec = RefStruct.particle_index_ref(frame_set_filter); 
        snip_data.nc_master_vec = RefStruct.master_nucleusID_ref(frame_set_filter);    

        % specify name      
        snip_name = ['snip_data_F' sprintf('%03d',samplingInfo.Frame) '_S' sprintf('%02d',samplingInfo.SetID)]; 
    %     if write_snip_flag            
    %         blank = struct;
    %         save([dataPath 'snip_data.mat'],'blank','-v7.3')    
    %         write_snip_flag = false;
    %     end
        snip_file = matfile([liveProject.dataPath 'snip_data.mat'],'Writable',true);    
        snip_file.(snip_name)= snip_data;        
    %     clear snip_file;    

        % update waitbar
        send(D, i);
    end


    % disp('saving qc frames...')
    % % save qc data
    % tic
    % particle_index = unique([spot_struct_protein.particleID]);
    % particle_index = particle_index(~isnan(particle_index));
    % qc_particles = randsample(particle_index,min([100,numel(particle_index)]),false);
    % particle_index_full = [];
    % particle_frames_full = [];
    % for i = 1:numel(qc_structure)
    %     qc_mat = qc_structure(i).qc_mat;
    %     for  j = 1:numel(qc_mat)
    %         qc_spot = qc_mat(j);
    %         if ~isfield(qc_spot,'ParticleID')
    %             continue
    %         end
    %         ParticleID = qc_spot.particleID;
    %         if isempty(ParticleID) || ~ismember(ParticleID,qc_particles)
    %             continue
    %         end        
    %         samplingInfo.Frame = qc_spot.frame;      
    %         particle_index_full = [particle_index_full ParticleID];
    %         particle_frames_full = [particle_frames_full samplingInfo.Frame];        
    %         save_name = [snipPath 'pt' num2str(1e4*ParticleID) '_frame' sprintf('%03d',samplingInfo.Frame) '.mat'];
    %         save(save_name,'qc_spot');
    %     end
    % end
    % [particle_index_full, si] = sort(particle_index_full);
    % particle_frames_full = particle_frames_full(si);
    % 
    % qc_ref_struct.particle_frames_full = particle_frames_full;
    % qc_ref_struct.particle_index_full = particle_index_full;
    % toc
    % save updated nucleus structure
    disp('saving nucleus structure...')

    % save([dataPath 'qc_ref_struct.mat'],'qc_ref_struct')
    save([liveProject.dataPath 'spot_struct_protein.mat'],'spot_struct_protein','-v7.3') 
    save([liveProject.dataPath 'proteinSamplingInfo.mat'],'proteinSamplingInfo');


%% %%%%%%%%% waitbar function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   function nUpdateWaitbar(~)
%       waitbar(p/N, h);
%       p = p + 1;
%   end
end