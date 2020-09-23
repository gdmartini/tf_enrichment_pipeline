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
    
    %% %%%%%%%%%%%%%%%%%%%%%%% Initialize enrichment-related fields  %%%%%%%%%%

    [spot_struct_protein, ~] = initializeProteinFields(spot_struct_protein, use3DSpotInfo);

    snip_fields = {'spot_protein_snips', 'edge_control_protein_snips',...
        'spot_mcp_snips','edge_control_mcp_snips'};

    % get list of protein-specific fields that were added
    spot_fields = fieldnames(spot_struct_protein);
    NewFields = spot_fields(~ismember(spot_fields,fieldnames(nucleus_struct)));
    
    %% %%%%%%%%%%%%%%%%%%%%%%% Generate indexing vectors  %%%%%%%%%%%%%%%%%%%%%
    [RefStruct, SetFrameArray, SamplingResults] = generateReferenceVectors(spot_struct_protein,refPath,use3DSpotInfo,ignoreQC,NewFields);

    qc_structure = struct;   

    %% %%%%%%%%%%%%%%%%%%%%%%% Nucleus segmentation  %%%%%%%%%%%%%%%%%%%%%%%%%%

    % first check to see if segmentation files exist
    [segmentNuclei, segmentIndices] = ...
              handleSegmentationOptions(RefStruct,segmentNuclei);

    if segmentNuclei
        disp('segmenting nuclei...')    
        nuclearSegmentation(liveProject, RefStruct, segmentIndices, NumWorkers);      
    end

    %% %%%%%%%%%%%%%%%% Local Protein Sampling: Loop 1 %%%%%%%%%%%%%%%%%%%%
    % Generate positions for serialized and edge-controlled computational
    % control spots. This can be done without loading the actual image
    % stacks so it is quite fast
        
    % convert structure to cell array for conveneince
    pool = gcp('nocreate');
    if isempty(pool)
      parpool
    end        
    
    
        D = parallel.pool.DataQueue;    
    afterEach(D, @nUpdateWaitbar);
    
%     N = size(set_frame_array,1);
%     p = 1;
% 
    h = waitbar(0,'Generating control spots...');
    
    for i = 1:size(SetFrameArray,1)  
        waitbar(i/size(SetFrameArray,1),h)       
        samplingInfo = struct; % use this structure to keep track of key tracking-related info

        % read basic info from set_frame array
        samplingInfo.SetID = SetFrameArray(i,1);
        samplingInfo.Frame = SetFrameArray(i,2);  

        samplingInfo = getSamplingInfo(samplingInfo,liveProject,proteinSamplingInfo,RefStruct);
        
        % End identical chunk
        
        % perform QC and generate lookup table of inter-nucleus distances
        samplingInfo = performNucleusQC(samplingInfo);

        j_pass = 1; % counter to track absolute position in iteration
        for j = samplingInfo.frame_set_indices        
            % get indexing info         
            samplingInfo.spotIndex = RefStruct.particle_index_ref(j);
            samplingInfo.spotSubIndex = RefStruct.particle_subindex_ref(j);

            x_spot = RefStruct.spot_x_ref(j);
            x_index = min([samplingInfo.xDim max([1 round(x_spot)])]);
            y_spot = RefStruct.spot_y_ref(j);
            y_index = min([samplingInfo.yDim max([1 round(y_spot)])]);

            % extract mask 
            nucleus_mask_id = samplingInfo.nc_label_frame(y_index,x_index);            
            
            % if spot does not fall within boundaries of a nucleus, flag it
            % and skip
            if ~nucleus_mask_id

                SamplingResults(i).edge_qc_flag_vec(j_pass) = -1;            
                SamplingResults(i).serial_qc_flag_vec(j_pass) = -1;

                continue
            end 
            
            % create mask
            nucleus_mask = samplingInfo.nc_label_frame == nucleus_mask_id; 
                 
            %% %%%% Find edge control sample location %%%%%%%%%%%%%%%%%%%%%
            [SamplingResults(i).edge_null_nc_vec(j_pass), ...
              SamplingResults(i).edge_null_x_vec(j_pass),...
              SamplingResults(i).edge_null_y_vec(j_pass), ...
              SamplingResults(i).edge_qc_flag_vec(j_pass)] = ...
              ...
              findEdgeControlWrapper(...
                samplingInfo,nucleus_mask,x_index,y_index,nucleus_mask_id);
            
            %% %%%% Find serialized control ocation %%%%%%%%%%%%%%%%%%%%%%%
            % This is the bit that cannot be easily parallelized
            [spot_struct_protein(samplingInfo.spotIndex).serial_null_edge_dist_vec(samplingInfo.spotSubIndex),...
             spot_struct_protein(samplingInfo.spotIndex).serial_null_x_vec(samplingInfo.spotSubIndex),...
             spot_struct_protein(samplingInfo.spotIndex).serial_null_y_vec(samplingInfo.spotSubIndex),...
             spot_struct_protein(samplingInfo.spotIndex).serial_qc_flag_vec(samplingInfo.spotSubIndex)]...
             ...
              = drawSerializedControlSpot(...
                  samplingInfo,nucleus_mask,spot_struct_protein(samplingInfo.spotIndex).frames, ...
                  spot_struct_protein(samplingInfo.spotIndex).serial_null_x_vec, ...
                  spot_struct_protein(samplingInfo.spotIndex).serial_null_y_vec);
                
             % pass spot structure values to sample structure
             SamplingResults(i).serial_null_edge_dist_vec(j_pass) = spot_struct_protein(samplingInfo.spotIndex).serial_null_edge_dist_vec(samplingInfo.spotSubIndex);
             SamplingResults(i).serial_null_x_vec(j_pass) = spot_struct_protein(samplingInfo.spotIndex).serial_null_x_vec(samplingInfo.spotSubIndex);
             SamplingResults(i).serial_null_y_vec(j_pass) = spot_struct_protein(samplingInfo.spotIndex).serial_null_y_vec(samplingInfo.spotSubIndex);
             SamplingResults(i).serial_qc_flag_vec(j_pass) = spot_struct_protein(samplingInfo.spotIndex).serial_qc_flag_vec(samplingInfo.spotSubIndex);
             
             % increment
             j_pass = j_pass + 1;
        end      
    end
    
    %% %%%%%%%%%%%%%%%% Local Protein Sampling: Loop 2 %%%%%%%%%%%%%%%%%%%%
    
    
    for i = 1:size(SetFrameArray,1)
        
        % generate structure to keep track of sampling info
        samplingInfo = struct; 
        
        samplingInfo.SetID = SetFrameArray(i,1);
        samplingInfo.Frame = SetFrameArray(i,2);  

        samplingInfo = getSamplingInfo(samplingInfo,liveProject,proteinSamplingInfo,RefStruct);             

        % load stacks            
        samplingInfo.protein_stack = imreadStack(samplingInfo.proteinPath);
        if size(samplingInfo.protein_stack,3) == samplingInfo.zDim+2
          samplingInfo.protein_stack = samplingInfo.protein_stack(:,:,2:end); 
        else
          error('Unrecognized form of z-padding')
        end    
       
        samplingInfo.mcp_stack = imreadStack2(samplingInfo.mcpPath);
        if size(samplingInfo.mcp_stack,3) == samplingInfo.zDim+2
          samplingInfo.mcp_stack = samplingInfo.mcp_stack(:,:,2:end);
        else
          error('Unrecognized form of z-padding')
        end

        % perform QC and generate lookup table of inter-nucleus distances
        samplingInfo = performNucleusQC(samplingInfo);           

        % initialize temporary arrays to store snip info 
        temp_snip_struct = struct;
        for j = 1:length(snip_fields)
            temp_snip_struct.(snip_fields{j})  = NaN(2*samplingInfo.snippet_size+1,2*samplingInfo.snippet_size+1,length(samplingInfo.frame_set_indices));
        end    
        
        % iterate through spots
        qc_mat = struct;
        
        j_pass = 1; % counter to track absolute position in iteration
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
            nucleus_mask_id = samplingInfo.nc_label_frame(y_index,x_index);                                                

            %% %%%%%%%%%%%%%%%%%%% Sample protein levels %%%%%%%%%%%%%%%%%%%%%%                                     

            % make sure size is reasonable and that spot is inside nucleus
            if ~nucleus_mask_id

                SamplingResults(i).edge_qc_flag_vec(j_pass) = -1;            
                SamplingResults(i).serial_qc_flag_vec(j_pass) = -1;
                
                nucleus_mask_3D_dummy = true(size(samplingInfo.protein_stack)); 
                
                % still sample protein near locus in this case
                SamplingResults(i).spot_protein_vec(j_pass) = ...
                  sample_protein_3D(samplingInfo,samplingInfo.protein_stack,...
                  x_spot,y_spot,z_spot,samplingInfo.xy_sigma,samplingInfo.z_sigma,nucleus_mask_3D_dummy);

                SamplingResults(i).spot_mcp_vec(j_pass) = ...
                  sample_protein_3D(samplingInfo,samplingInfo.mcp_stack,...
                  x_spot,y_spot,z_spot,samplingInfo.xy_sigma,samplingInfo.z_sigma,nucleus_mask_3D_dummy);

                continue
            end 
            
            % create mask
            nucleus_mask_3D = repmat(samplingInfo.nc_label_frame == nucleus_mask_id,size(samplingInfo.protein_stack,3)); 
            
            % sample protein near locus
            SamplingResults(i).spot_protein_vec(j_pass) = ...
              sample_protein_3D(samplingInfo,samplingInfo.protein_stack,...
              x_spot,y_spot,z_spot,samplingInfo.xy_sigma,samplingInfo.z_sigma,nucleus_mask_3D);

            SamplingResults(i).spot_mcp_vec(j_pass) = ...
              sample_protein_3D(samplingInfo,samplingInfo.mcp_stack,...
              x_spot,y_spot,z_spot,samplingInfo.xy_sigma,samplingInfo.z_sigma,nucleus_mask_3D);
            
            % sample snippets        
            temp_snip_struct.spot_protein_snips(:,:,j_pass) = sample_snip_3D(x_spot,y_spot,z_spot,samplingInfo,samplingInfo.protein_stack,nucleus_mask_3D);
            temp_snip_struct.spot_mcp_snips(:,:,j_pass) = sample_snip_3D(x_spot,y_spot,z_spot,samplingInfo,samplingInfo.mcp_stack,nucleus_mask_3D); 

            % Take average across all pixels within 1.5um of nuclues center           
            SamplingResults(i).nucleus_protein_vec(j_pass) = sample_protein_3D(...
              samplingInfo,samplingInfo.protein_stack,x_nucleus,y_nucleus,z_spot,samplingInfo.xy_sigma_nuclear,samplingInfo.z_sigma,nucleus_mask_3D);        

            %% %%%%%%%%%%%% Draw edge control spot %%%%%%%%%%%%%%%%%%%%%%%%%%%%                           

            % Draw control samples (as appropriate)    
            if SamplingResults(i).edge_qc_flag_vec(j_pass) > 0                
                edge_control_x = SamplingResults(i).edge_null_x_vec(j_pass);
                edge_control_y = SamplingResults(i).edge_null_y_vec(j_pass);     
                
                % get mask
                if SamplingResults(i).edge_qc_flag_vec(j_pass) == 1
                    nucleus_mask_3D_edge = nucleus_mask_3D;
                elseif SamplingResults(i).edge_qc_flag_vec(j_pass) == 2
                    % extract mask 
                    nn_nucleus_mask_id = samplingInfo.nc_label_frame(edge_control_y,edge_control_x); 
                    nucleus_mask_3D_edge = repmat(samplingInfo.nc_label_frame == nn_nucleus_mask_id,size(samplingInfo.protein_stack,3)); 
                end
                
                SamplingResults(i).edge_null_protein_vec(j_pass) = sample_protein_3D(...
                  samplingInfo,samplingInfo.protein_stack,edge_control_x,edge_control_y,z_spot,...
                  samplingInfo.xy_sigma,samplingInfo.z_sigma,nucleus_mask_3D_edge);

                SamplingResults(i).edge_null_mcp_vec(j_pass) = sample_protein_3D(...              
                  samplingInfo,samplingInfo.mcp_stack,edge_control_x,edge_control_y,z_spot,...
                  samplingInfo.xy_sigma,samplingInfo.z_sigma,nucleus_mask_3D_edge);

                % draw snips    
                temp_snip_struct.edge_control_protein_snips(:,:,j_pass) = sample_snip_3D(edge_control_x,edge_control_y,z_spot,samplingInfo,samplingInfo.protein_stack,nucleus_mask_3D_edge);
                temp_snip_struct.edge_control_mcp_snips(:,:,j_pass) = sample_snip_3D(edge_control_x,edge_control_y,z_spot,samplingInfo,samplingInfo.mcp_stack,nucleus_mask_3D_edge);
            end                  

            %% %%%%%%%%%%%% Draw serialized control spot %%%%%%%%%%%%%%%%%%%%%%%                                                                        
            
            % sample protein 
            serial_control_x = SamplingResults(i).serial_null_x_vec(j_pass);
            serial_control_y = SamplingResults(i).serial_null_y_vec(j_pass);

            SamplingResults(i).serial_null_protein_vec(j_pass) = sample_protein_3D(...
                  samplingInfo,samplingInfo.protein_stack,serial_control_x,serial_control_y,z_spot,samplingInfo.xy_sigma,samplingInfo.z_sigma,nucleus_mask_3D);
            SamplingResults(i).serial_null_mcp_vec(j_pass) = sample_protein_3D(...
                  samplingInfo,samplingInfo.mcp_stack,serial_control_x,serial_control_y,z_spot,samplingInfo.xy_sigma,samplingInfo.z_sigma,nucleus_mask_3D);


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
        snip_data.particle_id_vec = RefStruct.particleID_ref(FrameSetFilter);
    
        % indexing vectors    
        snip_data.nc_sub_index_vec = RefStruct.particle_subindex_ref(FrameSetFilter);
        snip_data.nc_lin_index_vec = RefStruct.particle_index_ref(FrameSetFilter); 
        snip_data.nc_master_vec = RefStruct.master_nucleusID_ref(FrameSetFilter);    

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