function [selectControlSpots, SamplingResults] = handleControlSelectionOptions(...
          proteinSamplingInfoNew,SamplingResults,selectControlSpots,liveProject)
        
    load([liveProject.dataPath 'spot_struct_protein.mat'],'spot_struct_protein','-v7.3') 
    load([liveProject.dataPath 'proteinSamplingInfo.mat'],'proteinSamplingInfo');


% if all(spot_frame_vec) && segmentNuclei   
%     warning('previous segmentation results found')
%     y = 1;
%     n = 0;
%     overwrite = input('overwrite segmentation results? (y/n)');
%     segmentNuclei = overwrite;
%     
% elseif ~all(nc_frame_vec) && ~segmentNuclei   
%     warning('some or all frames missing nucleus segmentation data. Segmenting missing frames only')
%     segmentNuclei = 1;   
%     segmentIndices = find(~nc_frame_vec);
% end