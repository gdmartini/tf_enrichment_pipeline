function spot_struct_protein = drawSerializedControlSpot(...
                                  spot_struct_protein,samplingInfo,...
                                x_index,y_index,spotIndex,spotSubIndex,nucleus_mask)                                  

  % spot and nucleus distance info
  spot_edge_dist = samplingInfo.nc_dist_frame(y_index,x_index);        
  nc_edge_dist_vec = samplingInfo.nc_dist_frame(nucleus_mask);
  
  spot_struct_protein(spotIndex).spot_edge_dist_vec(spotSubIndex) = spot_edge_dist;        
  spot_sep_vec = samplingInfo.spot_dist_frame(nucleus_mask); 
  x_pos_vec = samplingInfo.x_ref(nucleus_mask);
  y_pos_vec = samplingInfo.y_ref(nucleus_mask);         
  
  % time series info
  frame_vec = spot_struct_protein(spotIndex).frames; 
  currentFrame = double(frame_vec(spotSubIndex));
  serial_null_x = spot_struct_protein(spotIndex).serial_null_x_vec;
  serial_null_y = spot_struct_protein(spotIndex).serial_null_y_vec; 

  % generate sampling vector      
  sample_index_vec = find(spot_sep_vec >= samplingInfo.minSampleSep & ...
              nc_edge_dist_vec >= samplingInfo.minEdgeSep);
  
  % if this is the first sample for this spot, just find random
  % control snip. This will "seed" subsequent samples
  if ~isempty(sample_index_vec)
      if all(isnan(serial_null_x))          
    
          % Take a random sample filter for regions far enough away from locus    
      
          if length(sample_index_vec) > 1
            new_index = randsample(sample_index_vec,1);
          else
            new_index = sample_index_vec;
          end
          
          serial_control_x = x_pos_vec(new_index);
          serial_control_y = y_pos_vec(new_index);
          serial_edge_dist = nc_edge_dist_vec(new_index);                             
            
    % otherwise, draw snip based on previous location
    else
        prevIndex = (find(~isnan(serial_null_x),1,'last'));
        n_frames = currentFrame - double(frame_vec(prevIndex)); % used to adjust jump weights
        old_x = double(serial_null_x(prevIndex));
        old_y = double(serial_null_y(prevIndex));    

        % calculate distance from previous location     
        drControl = double(sqrt((old_x-x_pos_vec(sample_index_vec)).^2+(old_y-y_pos_vec(sample_index_vec)).^2));   

        % calculate weights
        wt_vec = exp(-.5*((drControl/double(sqrt(n_frames)*(samplingInfo.driftTol))).^2));      

        % draw sample
        if any(wt_vec>0)
            new_index = randsample(sample_index_vec,1,true,wt_vec);
            serial_control_x = x_pos_vec(new_index);
            serial_control_y = y_pos_vec(new_index);
            serial_edge_dist = nc_edge_dist_vec(new_index);  
            
        else
            new_index = randsample(sample_index_vec,1,true);
            serial_control_x = x_pos_vec(new_index);
            serial_control_y = y_pos_vec(new_index);
            serial_edge_dist = nc_edge_dist_vec(new_index);  
        end
      end
  else
      warning('Unable to draw serial control spot. Check nucleus segmentation, and "minEdgeSep" and "minSampSep" parameters')
  end
  
  % record info
  spot_struct_protein(spotIndex).serial_null_x_vec(spotSubIndex) = serial_control_x;
  spot_struct_protein(spotIndex).serial_null_y_vec(spotSubIndex) = serial_control_y;
  spot_struct_protein(spotIndex).serial_null_edge_dist_vec(spotSubIndex) = serial_edge_dist;