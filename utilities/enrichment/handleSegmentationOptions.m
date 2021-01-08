function [segmentNuclei, segmentIndices] = ...
          handleSegmentationOptions(refVecStruct,segmentNuclei)

set_frame_array = refVecStruct.set_frame_array;
refPath = refVecStruct.refPath;

segmentIndices = 1:size(set_frame_array,1);
spot_frame_vec = false(1,size(set_frame_array,1));
nc_frame_vec = false(1,size(set_frame_array,1));
for i = 1:size(set_frame_array,1)        
    setID_temp = set_frame_array(i,1);
    frame_temp = set_frame_array(i,2);  
    nc_ref_name = [refPath 'nc_ref_frame_set' sprintf('%02d',setID_temp) '_frame' sprintf('%03d',frame_temp) '.mat'];
    nc_frame_vec(i) = isfile(nc_ref_name);    
    spot_ref_name = [refPath 'spot_roi_frame_set' sprintf('%02d',setID_temp) '_frame' sprintf('%03d',frame_temp) '.mat'];    
    spot_frame_vec(i) = isfile(spot_ref_name);    
end    

if all(spot_frame_vec) && segmentNuclei   
    warning('previous segmentation results found')
    y = 1;
    n = 0;
    overwrite = input('overwrite segmentation results? (y/n)');
    segmentNuclei = overwrite;
    
elseif ~all(nc_frame_vec) && ~segmentNuclei   
    warning('some or all frames missing nucleus segmentation data. Segmenting missing frames only')
    segmentNuclei = 1;   
    segmentIndices = find(~nc_frame_vec);
end