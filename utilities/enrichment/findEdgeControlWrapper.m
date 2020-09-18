function [spot_edge_dist, edge_x, edge_y, qc_flag] = findEdgeControlWrapper(...
                                RefStruct,samplingInfo,nucleus_mask,x_index,y_index,...
                                spotIndex,spotSubIndex,j_pass)

  
  % Edge sampling 
  spot_edge_dist = samplingInfo.nc_dist_frame(y_index,x_index);        
  nc_edge_dist_vec = samplingInfo.nc_dist_frame(nucleus_mask);  
  spot_sep_vec = samplingInfo.spot_dist_frame(nucleus_mask);        

  % Now find control "spot" that is same distance from nucleus edge
  % as true spot
  [edge_x, edge_y, nucleus_id, qc_flag,~]...
      = find_control_sample(nc_edge_dist_vec, samplingInfo.x_ref, samplingInfo.y_ref, spot_sep_vec, spot_edge_dist,...
           j_pass, samplingInfo.minSampleSep, nucleus_mask,0);  

  % if initial attempt failed, try nearest neighbor nucleus
%   null_mask = spot_nc_mask;
  if qc_flag == 0

      % Find nearest neighbor nucleus
      r_vec = samplingInfo.r_dist_mat(:,j_pass);
      r_vec(j_pass) = Inf;
      [~, closestIndex] = min(r_vec);

      % get nn nucleus mask   
      x_spot_nn = RefStruct.spot_samplingInfo.x_ref(samplingInfo.frame_set_indices(closestIndex));
      y_spot_nn = RefStruct.spot_samplingInfo.y_ref(samplingInfo.frame_set_indices(closestIndex)); 

      nn_nucleus_mask = nc_ref_frame == RefStruct.master_nucleusID_ref(samplingInfo.frame_set_indices(closestIndex));
      
      if ~isnan(x_spot_nn)
          nan_flag = isnan(nn_nucleus_mask(y_spot_nn,x_spot_nn));
      end

      % make sure size is reasonable 
      if sum(nn_nucleus_mask(:)) >= sampleInfo.min_nucleus_area && sum(nn_nucleus_mask(:)) <= sampleInfo.max_nucleus_area && ~nan_flag                

          nn_edge_dist_vec = samplingInfo.nc_dist_frame(nn_nucleus_mask);
          nn_sep_vec = samplingInfo.spot_dist_frame(nn_nucleus_mask);

          [edge_x, edge_y, nucleus_id, qc_flag,~]...           
              = find_control_sample(nn_edge_dist_vec, samplingInfo.x_ref, samplingInfo.y_ref, nn_sep_vec, spot_edge_dist,...
                   closestIndex, samplingInfo.minSampleSep, nn_nucleus_mask,1);
                 
          if qc_flag == 1
              qc_flag = 2;
          end
      end
  end 
  