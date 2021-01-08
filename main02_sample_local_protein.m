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
%         snip_struct: data structure containing "snips" of nerhborhood
%         around sampled spots

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
    max_dist_nearest_neighbor_um = 2.5*max_nucleus_radius_um;
    segmentNuclei = 0;
    use3DSpotInfo = 0;
%     NumWorkers = 24;
    % PSF info for 3D sampling
    use_psf_fit_dims = false; % NL: currently no supported
    xy_sigma_um = 0.25;% um 
    xy_sigma_nuclear_um = 1.5;
    z_sigma_um = 0.6; % um
    ignoreQC = true;    
    NumWorkers = [];
    segmentationMethod = 1;
    overwriteSegmentation = false;

    %% %%%%%%%%%%%%%%%%%%%%%%% Check for optional inputs %%%%%%%%%%%%%%%%%%%%%%

    %options must be specified as name, value pairs. unpredictable errors will
    %occur, otherwise.
    for i = 1:2:(numel(varargin)-1)
        if i ~= numel(varargin)
            eval([varargin{i} '=varargin{i+1};']);
        end
    end

    parDefaultFlag = isempty(NumWorkers);
    
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
    proteinSamplingInfo.max_dist_nearest_neighbor_um = max_dist_nearest_neighbor_um;
    
    %% %%%%%%%%%%%%%%%%%%%%%%% Get project info %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [liveProject, ~, dataName, hasAPInfo, has3DSpotInfo, hasProteinInfo, hasNucleusProbFiles] = headerFunction(projectName);
    if ~hasProteinInfo
      warning('No input protein info associated with this project. Aborting protein sampling')
      return
    end
    use3DSpotInfo = use3DSpotInfo&&has3DSpotInfo;
    proteinSamplingInfo.use3DSpotInfo = use3DSpotInfo;

    %% %%%%%%%%%%%%%%%%%%%%%%% Load data and clean trace data %%%%%%%%%%%%%%%%%
    load(dataName,'spot_struct')

    if false % NL: Currently not working
        load([liveProject.dataPath '/psf_dims.mat'],'psf_dims')
    end

    % make paths 
    snipPathTemp = [liveProject.dataPath 'snip_fragments_temp\'];
    qcPath = [liveProject.dataPath 'qc_images\'];
    refPath = [liveProject.dataPath 'refFrames\'];

    % no make the directories
    mkdir(refPath)
    mkdir(qcPath)
    mkdir(snipPathTemp);    
    
    % remove old protein snip files (if any exist)
    snipCleanup(snipPathTemp)
    
    % remove frames where no particle was observed
    spot_struct_protein = truncateParticleFields(spot_struct,use3DSpotInfo,hasAPInfo);  
    
    %% %%%%%%%%%%%%%%%%%%%%%%% Initialize enrichment-related fields  %%%%%%%%%%

    [spot_struct_protein, ~] = initializeProteinFields(spot_struct_protein, use3DSpotInfo);

    snip_fields = {'spot_protein_snips', 'edge_control_protein_snips',...
        'spot_mcp_snips','edge_control_mcp_snips'};

    % get list of protein-specific fields that were added
    spot_fields = fieldnames(spot_struct_protein);
    NewFields = spot_fields(~ismember(spot_fields,fieldnames(spot_struct)));
    
    %% %%%%%%%%%%%%%%%%%%%%%%% Generate indexing vectors  %%%%%%%%%%%%%%%%%%%%%
    [RefStruct, SetFrameArray, SamplingResults] = generateReferenceVectors(...
      spot_struct_protein,refPath,use3DSpotInfo,ignoreQC,NewFields);

    qc_structure = struct;   

    %% %%%%%%%%%%%%%%%%%%%%%%% Nucleus segmentation  %%%%%%%%%%%%%%%%%%%%%%%%%%

    % first check to see if segmentation files exist
    [segmentNuclei, segmentIndices] = ...
              handleSegmentationOptions(RefStruct,segmentNuclei);

    
    if segmentNuclei
%         disp('segmenting nuclei...')   
        if ~parDefaultFlag
          nuclearSegmentation(liveProject, RefStruct, segmentIndices, spot_struct, NumWorkers, segmentationMethod, hasNucleusProbFiles);      
        else
          nuclearSegmentation(liveProject, RefStruct, segmentIndices, spot_struct, [], segmentationMethod, hasNucleusProbFiles);      
        end
    end

    %% %%%%%%%%%%%%%%%% Local Protein Sampling: Loop 1 %%%%%%%%%%%%%%%%%%%%
    % Generate positions for serialized and edge-controlled computational
    % control spots. This can be done without loading the actual image
    % stacks so it is quite fast
             
    h = waitbar(0,'Generating control spots...');        
    NIter = size(SetFrameArray,1);
    for i = 1:NIter%size(SetFrameArray,1)  
        waitbar(i/NIter,h)                       
        
        NewSetFlag = 1;
        if i == 1
          samplingInfo = struct; % use this structure to keep track of key tracking-related info
        elseif i > 1 && SetFrameArray(i,1)~=SetFrameArray(i-1,1)     
          samplingInfo = struct;
        else
          NewSetFlag = 0;
        end
        
        % read basic info from set_frame array
        samplingInfo.SetID = SetFrameArray(i,1);
        samplingInfo.Frame = SetFrameArray(i,2); 
        
        samplingInfo = getSamplingInfo(samplingInfo,liveProject,proteinSamplingInfo,RefStruct,NewSetFlag);                
        
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
                
                j_pass = j_pass + 1;
                continue
            end 
            
            % create mask
            nucleus_mask = samplingInfo.nc_label_frame == nucleus_mask_id; 
                 
            %% %%%% Find edge control sample location %%%%%%%%%%%%%%%%%%%%%
            [SamplingResults(i).spot_edge_dist_vec(j_pass), ...
              SamplingResults(i).edge_null_x_vec(j_pass),...
              SamplingResults(i).edge_null_y_vec(j_pass), ...
              SamplingResults(i).edge_qc_flag_vec(j_pass)] = ...
              ...
              findEdgeControlWrapper(...
                samplingInfo,nucleus_mask,x_index,y_index,nucleus_mask_id);
            
            %% %%%% Find serialized control location %%%%%%%%%%%%%%%%%%%%%%%
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
    delete(h);
    %% %%%%%%%%%%%%%%%% Local Protein Sampling: Loop 2 %%%%%%%%%%%%%%%%%%%%
    h = waitbar(0,'Sampling protein...'); 
    
    pool = gcp('nocreate');
    if isempty(pool)
      if parDefaultFlag
        parpool;
      else
        parpool(NumWorkers);
      end  
    end        
       
    D = parallel.pool.DataQueue;    
    afterEach(D, @nUpdateWaitbar);
    
    N = size(SetFrameArray,1);
    p = 1;
    
    tic
    for i = 1:NIter%size(SetFrameArray,1)
        % generate structure to keep track of sampling info
        samplingInfo = struct; 
        
        samplingInfo.SetID = SetFrameArray(i,1);
        samplingInfo.Frame = SetFrameArray(i,2);  

        samplingInfo = getSamplingInfo(samplingInfo,liveProject,proteinSamplingInfo,RefStruct,1);             
  
        % load stacks            
        samplingInfo.protein_stack = imreadStack(samplingInfo.proteinPath);
        % find and remove padding slices
        padding_flags_protein = reshape(max(max(samplingInfo.protein_stack,[],1),[],2),[],1)==0;        
        samplingInfo.protein_stack = samplingInfo.protein_stack(:,:,~padding_flags_protein); 
               
        samplingInfo.mcp_stack = imreadStack(samplingInfo.mcpPath);
        % find and remove padding slices
        padding_flags_mcp = reshape(max(max(samplingInfo.mcp_stack,[],1),[],2),[],1)==0;        
        samplingInfo.mcp_stack = samplingInfo.mcp_stack(:,:,~padding_flags_mcp); 
        if any(padding_flags_mcp~=padding_flags_protein)
            error('inconsistent padding between MCP and Protein channels')
        end
        
        % determine z offset
        zOffset = length(padding_flags_mcp)-find(~padding_flags_mcp,1,'last');
        
        % reset zDim parameter to ensure consistency
        if samplingInfo.zDim ~= sum(padding_flags_mcp)
            samplingInfo.zDim = sum(~padding_flags_mcp);
            [samplingInfo.x_ref,samplingInfo.y_ref,samplingInfo.z_ref] = meshgrid(1:samplingInfo.xDim,1:samplingInfo.yDim,1:samplingInfo.zDim);
        end
        % perform QC and generate lookup table of inter-nucleus distances
        samplingInfo = performNucleusQC(samplingInfo);           

        % initialize temporary arrays to store snip info 
        TempSnipStruct = struct;
        for j = 1:length(snip_fields)
            TempSnipStruct.(snip_fields{j})  = NaN(2*samplingInfo.snippet_size+1,2*samplingInfo.snippet_size+1,length(samplingInfo.frame_set_indices));
        end    
        
        % iterate through spots
        qc_mat = struct;
        
        j_pass = 1; % counter to track absolute position in iteration
        for j = samplingInfo.frame_set_indices 
            samplingSubInfo = struct;
            
            % get indexing info (is this used?)         
            samplingSubInfo.spotIndex = RefStruct.particle_index_ref(j);
            samplingSubInfo.spotSubIndex = RefStruct.particle_subindex_ref(j);

            % get location info
            samplingSubInfo.x_spot_full = RefStruct.spot_x_ref(j);
            samplingSubInfo.x_index = min([samplingInfo.xDim max([1 round(samplingSubInfo.x_spot_full)])]);
            samplingSubInfo.y_spot_full = RefStruct.spot_y_ref(j);
            samplingSubInfo.y_index = min([samplingInfo.yDim max([1 round(samplingSubInfo.y_spot_full)])]);
            samplingSubInfo.z_spot_full = RefStruct.spot_z_ref(j)-zOffset;%1.0; % adjust for z padding
            samplingSubInfo.z_index = min([samplingInfo.zDim max([1 round(samplingSubInfo.z_spot_full)])]);
            
            samplingSubInfo.z_vol_dim = ceil(2*samplingInfo.z_sigma); % size of sub-stack to sample
            samplingSubInfo.z_range3 = max(1,samplingSubInfo.z_index-samplingSubInfo.z_vol_dim):min(samplingInfo.zDim,samplingSubInfo.z_index+samplingSubInfo.z_vol_dim);
            samplingSubInfo.z_range3_full = 1:samplingInfo.zDim;
            
            % extract mask 
            nucleus_mask_id = samplingInfo.nc_label_frame(samplingSubInfo.y_index,samplingSubInfo.x_index);                                                

            %% %%%%%%%%%%%%%%%%%%% Sample protein levels %%%%%%%%%%%%%%%%%%%%%%                                     

            % make sure size is reasonable and that spot is inside nucleus
            if ~nucleus_mask_id

                SamplingResults(i).edge_qc_flag_vec(j_pass) = -1;            
                SamplingResults(i).serial_qc_flag_vec(j_pass) = -1;
                
                nucleus_mask_3D_dummy = true(size(samplingInfo.protein_stack)); 
                
                % still sample protein near locus in this case                
                SamplingResults(i).spot_protein_vec(j_pass) = ...
                  sample_protein_3D(samplingInfo,samplingInfo.protein_stack,...
                  samplingSubInfo.x_spot_full,samplingSubInfo.y_spot_full,...
                  samplingSubInfo.z_spot_full,samplingInfo.xy_sigma,samplingInfo.z_sigma,nucleus_mask_3D_dummy);
                
                SamplingResults(i).spot_mcp_vec(j_pass) = ...
                  sample_protein_3D(samplingInfo,samplingInfo.mcp_stack,...
                  samplingSubInfo.x_spot_full,samplingSubInfo.y_spot_full,...
                  samplingSubInfo.z_spot_full,samplingInfo.xy_sigma,samplingInfo.z_sigma,nucleus_mask_3D_dummy);
                
                % increment
                j_pass = j_pass + 1;
                
                continue
            end 
            
            % this generates smaller stack to draw samples from
            samplingSubInfo = createSubStack(samplingInfo,samplingSubInfo,nucleus_mask_id);            
            
            % sample protein near locus                     
            SamplingResults(i).spot_protein_vec(j_pass) = ...
              sample_protein_3D(samplingSubInfo,samplingSubInfo.protein_stack,...
              samplingSubInfo.x_spot,samplingSubInfo.y_spot,samplingSubInfo.z_spot,samplingInfo.xy_sigma,...
              samplingInfo.z_sigma,samplingSubInfo.nucleus_mask_3D);
                    
            SamplingResults(i).spot_mcp_vec(j_pass) = ...              
              sample_protein_3D(samplingSubInfo,samplingSubInfo.mcp_stack,...
              samplingSubInfo.x_spot,samplingSubInfo.y_spot,samplingSubInfo.z_spot,samplingInfo.xy_sigma,...
              samplingInfo.z_sigma,samplingSubInfo.nucleus_mask_3D);
            
            % sample snippets        
            TempSnipStruct.spot_protein_snips(:,:,j_pass) = ...
              sample_snip_3D(samplingSubInfo.x_spot,samplingSubInfo.y_spot,samplingSubInfo.z_spot,...
              samplingSubInfo,samplingSubInfo.protein_stack,samplingSubInfo.nucleus_mask_3D);
            
            TempSnipStruct.spot_mcp_snips(:,:,j_pass) = ...
              sample_snip_3D(samplingSubInfo.x_spot,samplingSubInfo.y_spot,samplingSubInfo.z_spot,...
              samplingSubInfo,samplingSubInfo.mcp_stack,samplingSubInfo.nucleus_mask_3D);

            % Take average across all pixels within 1.5um of nuclues center           
            SamplingResults(i).nuclear_protein_vec(j_pass) = sample_protein_3D(...
              samplingSubInfo,samplingSubInfo.mcp_stack,...
              samplingSubInfo.x_nucleus,samplingSubInfo.y_nucleus,samplingSubInfo.z_spot,samplingInfo.xy_sigma_nuclear,...
              samplingInfo.z_sigma,samplingSubInfo.nucleus_mask_3D);

            %% %%%%%%%%%%%% Draw edge control spot %%%%%%%%%%%%%%%%%%%%%%%%%%%%                           

            % Draw control samples (as appropriate)    
            if SamplingResults(i).edge_qc_flag_vec(j_pass) > 0                
                edge_control_x = SamplingResults(i).edge_null_x_vec(j_pass);
                edge_control_y = SamplingResults(i).edge_null_y_vec(j_pass);     
                
                % get mask
                if SamplingResults(i).edge_qc_flag_vec(j_pass) == 1
                    samplingSubInfoNN = samplingSubInfo;
                elseif SamplingResults(i).edge_qc_flag_vec(j_pass) == 2
                    % extract mask 
                    nn_nucleus_mask_id = samplingInfo.nc_label_frame(edge_control_y,edge_control_x); 
                    samplingSubInfoNN = createSubStack(samplingInfo,samplingSubInfo,nn_nucleus_mask_id);                    
                end                
                
                edge_control_x = edge_control_x + samplingSubInfoNN.x_shift;
                edge_control_y = edge_control_y + samplingSubInfoNN.y_shift;
                
                SamplingResults(i).edge_null_protein_vec(j_pass) = ...
                  sample_protein_3D(samplingSubInfoNN,samplingSubInfoNN.protein_stack,...
                  edge_control_x,edge_control_y,samplingSubInfoNN.z_spot,samplingInfo.xy_sigma,...
                  samplingInfo.z_sigma,samplingSubInfoNN.nucleus_mask_3D);
               
                SamplingResults(i).edge_null_mcp_vec(j_pass) = ...
                  sample_protein_3D(samplingSubInfoNN,samplingSubInfoNN.mcp_stack,...
                  edge_control_x,edge_control_y,samplingSubInfoNN.z_spot,samplingInfo.xy_sigma,...
                  samplingInfo.z_sigma,samplingSubInfoNN.nucleus_mask_3D);

                % draw snips    
                TempSnipStruct.edge_control_protein_snips(:,:,j_pass) = ...
                  sample_snip_3D(edge_control_x,edge_control_y,samplingSubInfoNN.z_spot,...
                  samplingSubInfoNN,samplingSubInfoNN.protein_stack,samplingSubInfoNN.nucleus_mask_3D);
                
                TempSnipStruct.edge_control_mcp_snips(:,:,j_pass) = ...
                  sample_snip_3D(edge_control_x,edge_control_y,samplingSubInfoNN.z_spot,...
                  samplingSubInfoNN,samplingSubInfoNN.mcp_stack,samplingSubInfoNN.nucleus_mask_3D);
            end                  

            %% %%%%%%%%%%%% Draw serialized control spot %%%%%%%%%%%%%%%%%%%%%%%                                                                                
            if SamplingResults(i).serial_qc_flag_vec(j_pass) > 0      
                % sample protein 
                serial_control_x = SamplingResults(i).serial_null_x_vec(j_pass) + samplingSubInfo.x_shift;                
                serial_control_y = SamplingResults(i).serial_null_y_vec(j_pass) + samplingSubInfo.y_shift;

                SamplingResults(i).serial_null_protein_vec(j_pass) = ...
                      sample_protein_3D(samplingSubInfo,samplingSubInfo.protein_stack,...
                      serial_control_x,serial_control_y,samplingSubInfo.z_spot,samplingInfo.xy_sigma,...
                      samplingInfo.z_sigma,samplingSubInfo.nucleus_mask_3D);

                SamplingResults(i).serial_null_mcp_vec(j_pass) = ...
                      sample_protein_3D(samplingSubInfo,samplingSubInfo.mcp_stack,...
                      serial_control_x,serial_control_y,samplingSubInfo.z_spot,samplingInfo.xy_sigma,...
                      samplingInfo.z_sigma,samplingSubInfo.nucleus_mask_3D);
            end

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
        saveSnipData(samplingInfo,TempSnipStruct,RefStruct,liveProject)    

        % update waitbar
        send(D, i);
    end
    toc
    delete(h);
    % map sampled fields back to spot structure    
    spot_struct_protein = mapToSpotStructure(spot_struct_protein,SamplingResults,RefStruct,SetFrameArray);
    
    % assemble snip fragments into a single mat file
    assembleSnips(liveProject);            
    
    % qc_ref_struct.particle_frames_full = particle_frames_full;
    % qc_ref_struct.particle_index_full = particle_index_full;
    % toc
    
    % run some basic data quality checks
    samplingChecks = performBasicSampleQC(spot_struct_protein);
    
    % save updated nucleus structure
    disp('saving nucleus structure...')

    % save([dataPath 'qc_ref_struct.mat'],'qc_ref_struct')
    save([liveProject.dataPath 'spot_struct_protein.mat'],'spot_struct_protein','-v7.3') 
    save([liveProject.dataPath 'proteinSamplingInfo.mat'],'proteinSamplingInfo');
    save([liveProject.dataPath 'samplingChecks.mat'],'samplingChecks');


%% %%%%%%%%% waitbar function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function nUpdateWaitbar(~)
      waitbar(p/N, h);
      p = p + 1;
  end
end