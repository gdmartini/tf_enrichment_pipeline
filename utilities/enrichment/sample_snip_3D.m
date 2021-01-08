function snip2D = sample_snip_3D(x_samp,y_samp,z_samp,samplingInfo,data_stack,nucleus_mask_3D)
    
    % mask out nuclear exterior (need to make this truly 3 dimensional)
    data_stack(~nucleus_mask_3D) = NaN;
    
    % generate rounded pos indices
    x_round = round(x_samp);
    y_round = round(y_samp);
    z_round = round(z_samp);
    
    % calculate vol dimensions
    xy_vol_dim = samplingInfo.snippet_size;
    z_vol_dim = ceil(2*samplingInfo.z_sigma);
    
    % calculate dimensions
    xDim = size(samplingInfo.z_ref,2);
    yDim = size(samplingInfo.z_ref,1);
    zDim = size(samplingInfo.z_ref,3);
   
    % volume protein sampling 
    x_range3 = max(1,x_round-xy_vol_dim):min(xDim,x_round+xy_vol_dim);
    y_range3 = max(1,y_round-xy_vol_dim):min(yDim,y_round+xy_vol_dim);
    z_range3 = max(1,z_round-z_vol_dim):min(zDim,z_round+z_vol_dim);
    x_range3_full = x_round-xy_vol_dim:x_round+xy_vol_dim;
    y_range3_full = y_round-xy_vol_dim:y_round+xy_vol_dim;
    z_range3_full = z_round-z_vol_dim:z_round+z_vol_dim;

    % generate protein sample and ref boxes
    protein_samp_box = NaN(numel(y_range3_full),numel(x_range3_full),numel(z_range3_full));
    protein_samp_box(ismember(y_range3_full,y_range3),ismember(x_range3_full,x_range3),...
        ismember(z_range3_full,z_range3)) = data_stack(y_range3,x_range3,z_range3);
      
    % z ref    
    z_ref_box = NaN(numel(y_range3_full),numel(x_range3_full),numel(z_range3_full));
    z_ref_box(ismember(y_range3_full,y_range3),ismember(x_range3_full,x_range3),...
        ismember(z_range3_full,z_range3)) = samplingInfo.z_ref(y_range3,x_range3,z_range3);
      
    % generate weight box
    weight_box = exp(-.5*(((z_samp-z_ref_box)./samplingInfo.z_sigma).^2));
    weight_box(weight_box<exp(-1/2)) = 0;  
    
    % take weighted average
    snip2D = nansum(protein_samp_box.*weight_box,3) ./ nansum(~isnan(protein_samp_box).*weight_box,3);
    