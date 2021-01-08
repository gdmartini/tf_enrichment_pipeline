function [RefStruct, SetFrameArray, SamplingResults] = generateReferenceVectors(spot_struct_protein,refPath,use3DSpotInfo,ignoreQC,NewFields)

%% %%%%%%%%%%%%%%%%%% generate ref structure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize structure to store reference vectors 
RefStruct = struct;
RefStruct.refPath = refPath;

% All of these vectors have an element for each spot detection in the set
RefStruct.frame_ref = [spot_struct_protein.frames];
RefStruct.nc_x_ref = [spot_struct_protein.xPosNucleus];
RefStruct.nc_y_ref = [spot_struct_protein.yPosNucleus];

if use3DSpotInfo
    RefStruct.spot_x_ref = [spot_struct_protein.xPosParticle3D]; 
    RefStruct.spot_y_ref = [spot_struct_protein.yPosParticle3D];
    RefStruct.spot_z_ref = [spot_struct_protein.zPosParticle3D];
else
  RefStruct.spot_x_ref = [spot_struct_protein.xPosParticle];
  RefStruct.spot_y_ref = [spot_struct_protein.yPosParticle];
  RefStruct.spot_z_ref = [spot_struct_protein.zPosParticle];
end
  
RefStruct.setID_ref = [];
RefStruct.master_nucleusID_ref = []; % list of nucleus indices
RefStruct.particle_index_ref = []; % list of particle indices within nucleus_struct
RefStruct.particle_subindex_ref = []; % list of frame indices within each particle
RefStruct.particleID_ref = [];
RefStruct.particle_qcFlag_ref = [];
for i = 1:numel(spot_struct_protein)
    ParticleID = spot_struct_protein(i).particleID;    
    RefStruct.particleID_ref = [RefStruct.particleID_ref repelem(ParticleID, numel(spot_struct_protein(i).frames))];
    RefStruct.particle_qcFlag_ref = [RefStruct.particle_qcFlag_ref spot_struct_protein(i).FrameQCFlags];
    RefStruct.setID_ref = [RefStruct.setID_ref repelem(spot_struct_protein(i).setID,numel(spot_struct_protein(i).frames))];
    RefStruct.master_nucleusID_ref = [RefStruct.master_nucleusID_ref repelem(round(spot_struct_protein(i).ncID*1e5),numel(spot_struct_protein(i).frames))];
    RefStruct.particle_index_ref = [RefStruct.particle_index_ref repelem(i,numel(spot_struct_protein(i).frames))]; 
    RefStruct.particle_subindex_ref = [RefStruct.particle_subindex_ref 1:numel(spot_struct_protein(i).frames)];
end

% option to override qc 
if ignoreQC 
    RefStruct.particle_qcFlag_ref = true(size(RefStruct.particleID_ref));
end

%%% Generate reference array for set-frame combos
RefStruct.set_frame_array = unique([RefStruct.setID_ref' RefStruct.frame_ref'],'row');
SetFrameArray = RefStruct.set_frame_array;

%% %%%%%%%%%%%%%%%%%% Make temp structure to store output %%%%%%%%%%%%%%%%%
SamplingResults = struct;
for i = 1:size(SetFrameArray,1)
  % read basic info from set_frame array
    SetID = SetFrameArray(i,1);
    Frame = SetFrameArray(i,2);  
    FrameSetFilter = RefStruct.setID_ref==SetID&RefStruct.frame_ref==Frame;
    for nf = 1:length(NewFields)
      SamplingResults(i).(NewFields{nf}) = NaN(1,sum(FrameSetFilter));
    end
end