function nuclearSegmentation(liveProject, RefStruct, segmentIndices, nucleus_struct, ...
        NumWorkers, segmentationMethod, hasNucleusProbFiles)  
  

  addpath(genpath('P:/Armando/LivemRNA/mRNADynamics\src')); 

  if NumWorkers ~=1 
      p = gcp('nocreate');
      if ~isempty(NumWorkers)    
        if isempty(p) || p.NumWorkers~=NumWorkers
          parpool(NumWorkers);
        end
      elseif isempty(p)
        parpool;
      end
  end
  
  
  
  %% %%%%%%%%%%%%%%%%% Make nucleus indexing vectors %%%%%%%%%%%%%%%%%%%%%%
  % construct a few indexing vectors using the full set of nucleus data
  % (spot or no spot)
  nucleus_set_ref = [];
  nucleus_frame_ref = [];
  nucleus_x_ref = [];
  nucleus_y_ref = [];
  for i = 1:length(nucleus_struct)
    nucleus_frame_ref = [nucleus_frame_ref nucleus_struct(i).frames];
    nucleus_x_ref = [nucleus_x_ref nucleus_struct(i).xPosNucleus];
    nucleus_y_ref = [nucleus_y_ref nucleus_struct(i).yPosNucleus];
    nucleus_set_ref = [nucleus_set_ref repelem(nucleus_struct(i).setID,length(nucleus_struct(i).frames))];
  end
  
  %% %%%%%%%%%%%%%%%%% set up parfor waitbar %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  D = parallel.pool.DataQueue;    
  afterEach(D, @nUpdateWaitbar);

  N = length(segmentIndices);
  p = 1;
    
  h = waitbar(0,'Segmenting nuclei...'); 

  % initialize arrays to store segmentation info
  nucleus_frame_cell = cell(1,length(segmentIndices));
  spot_frame_cell = cell(1,length(segmentIndices));

  parfor w = 1:length(segmentIndices)
      subInd = segmentIndices(w);
      currentSetID = RefStruct.set_frame_array(subInd,1);
      currentFrame = RefStruct.set_frame_array(subInd,2);  
      
      % get nucleus
      frame_set_filter_spot = RefStruct.setID_ref==currentSetID&RefStruct.frame_ref==currentFrame;
      frame_set_filter_nucleus = nucleus_set_ref==currentSetID&nucleus_frame_ref==currentFrame;
      
      % Get experiment info
      Prefix = liveProject.includedExperimentNames{currentSetID}; 
      currExperiment = liveProject.includedExperiments{currentSetID};                  
      proteinChannel = currExperiment.inputChannels;
      if length(proteinChannel) > 1
        error(['Problem with Prefix: ' Prefix '. This pipeline does not currently support multiple input channels'])
      end     
      
      % Generate reference vectors
      xDim = currExperiment.xDim;
      yDim = currExperiment.yDim;
      zDim = currExperiment.zDim; 
      
      % get nucleus positions
      nc_x_vec = round(nucleus_x_ref(frame_set_filter_nucleus));
      nc_x_vec(nc_x_vec<1) = 1;
      nc_x_vec(nc_x_vec>xDim) = xDim;
      nc_y_vec = round(nucleus_y_ref(frame_set_filter_nucleus)); 
      nc_y_vec(nc_y_vec<1) = 1;
      nc_y_vec(nc_y_vec>yDim) = yDim;
      
      % indexing vectors    
%       nc_master_vec = RefStruct.master_nucleusID_ref(frame_set_filter);  
      
      % get list of unique indices 
%       [nc_master_vec_u,ia,~] = unique(nc_master_vec,'stable');
      
      % unique nucleus vectors
%       nc_x_vec_u = round(nc_x_vec_temp(ia));
%       nc_y_vec_u = round(nc_y_vec_temp(ia));    
      
      % particle positions        
      spot_x_vec = round(RefStruct.spot_x_ref(frame_set_filter_spot));
      spot_y_vec = round(RefStruct.spot_y_ref(frame_set_filter_spot));                        

      % git pixel size info        
      PixelSize = currExperiment.pixelSize_nm / 1e3; % convert to microns
      nucleus_neighborhood_size = round(10 ./ PixelSize);  % determine size of neighborhood to use during nucleus segmentation      
      smoothing_kernel_size = round(1 ./ PixelSize); % size of gaussian smoothing kernel 

      % get protein channel       
      stackPath = [currExperiment.preFolder  Prefix '_' sprintf('%03d',currentFrame) '_ch0' num2str(proteinChannel) '.tif'];
      segment_stack = imreadStack2(stackPath, yDim, xDim, zDim+2);      
      segment_stack = segment_stack(:,:,2:end-1);            

      
      
      if segmentationMethod == 1 % default method
            nc_ref_frame = generateNuclearMask(segment_stack, nucleus_neighborhood_size, smoothing_kernel_size, xDim, yDim,...
                nc_y_vec, nc_x_vec, PixelSize);
      elseif segmentationMethod == 2 % Armando's method
            his = getHisMat(currExperiment);
            im = his(:, :, currentFrame);
            nc_ref_frame = kSnakeCircles(im,...
                PixelSize, 'fitEllipses', false, 'shouldWatershed', false);
      end
      
      displayFigures = false;
      
      if displayFigures 
          fig = figure;
          tiledlayout('flow');
          nexttile;
          imshow(im, []);
          nexttile;
          imshow(nc_ref_frame, []);
          waitforbuttonpress;
          close(fig);
      end
      
      % generate array indicating distance of each pixel from an active locus 
      nc_indices = sub2ind(size(nc_ref_frame),spot_y_vec,spot_x_vec);
      spot_dist_frame_temp = zeros(size(nc_ref_frame));
      spot_dist_frame_temp(nc_indices(~isnan(nc_indices))) = 1;
      spot_dist_frame_temp = bwdist(spot_dist_frame_temp);
      
      % store results
      nucleus_frame_cell{w} = nc_ref_frame;
      spot_frame_cell{w} = spot_dist_frame_temp;      
      
      % update waitbar
      send(D, w);
  end
  delete(h);
  disp('saving segmentation results...')
  % save arrays
  for w = 1:numel(segmentIndices)
      i = segmentIndices(w);
      currentSetID = RefStruct.set_frame_array(i,1);
      currentFrame = RefStruct.set_frame_array(i,2);  
      nc_ref_name = [RefStruct.refPath 'nc_ref_frame_set' sprintf('%02d',currentSetID) '_frame' sprintf('%03d',currentFrame) '.mat'];
      nc_ref_frame = nucleus_frame_cell{w};
      save(nc_ref_name,'nc_ref_frame');
      spot_ref_name = [RefStruct.refPath 'spot_roi_frame_set' sprintf('%02d',currentSetID) '_frame' sprintf('%03d',currentFrame) '.mat'];
      spot_dist_frame = spot_frame_cell{w};
      save(spot_ref_name,'spot_dist_frame');
  end
  disp('done.')  
  
  cleanUpmRNADynamics;
  
  %% %%%%%%%%% waitbar function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function nUpdateWaitbar(~)
      waitbar(p/N, h);
      p = p + 1;
  end



end


function nuclearMask = generateNuclearMask(nuclearImage, nucleus_neighborhood_size, smoothing_kernel_size,...
      xDim, yDim, nc_y_vec, nc_x_vec, PixelSize, varargin)

      protein_stack = nuclearImage;
      
      % generate protein gradient frame for segmentation
      protein_smooth = imgaussfilt(mean(protein_stack,3),round(smoothing_kernel_size/2));                
      protein_grad = imgradient(protein_smooth);   
      
      % flatten background         
      protein_bkg = imgaussfilt(protein_smooth, round(nucleus_neighborhood_size/2));
      protein_grad_norm = protein_grad ./ protein_bkg;        

      % try local thresholding
      n_x_local = floor(xDim / nucleus_neighborhood_size);
      n_y_local = floor(yDim / nucleus_neighborhood_size);
      x_dim_local = round(xDim / n_x_local);
      y_dim_local = round(yDim / n_y_local);
      protein_bin_clean = false(size(protein_grad_norm));
      for x = 1:n_x_local
          for y = 1:n_y_local
              x_start = (x-1)*x_dim_local+1;
              x_stop = min([x*x_dim_local,xDim]);
              y_start = (y-1)*y_dim_local+1;
              y_stop = min([y*y_dim_local,yDim]);
              pt_section = protein_grad_norm(y_start:y_stop,x_start:x_stop);
              thresh = multithresh(pt_section);
              section_bin = pt_section > thresh; 
              section_bin_clean = bwareaopen(section_bin,smoothing_kernel_size^2);                
              % record
              protein_bin_clean(y_start:y_stop,x_start:x_stop) = bwmorph(section_bin_clean,'hbreak');
          end
      end        
      
      % label regions
      seDim = max([1,round(3*.1/PixelSize)]);
      se = strel('disk',seDim);
      seBig = strel('disk',seDim+1);
      
      nc_frame_mask = logical(protein_bin_clean);%imerode(logical(protein_bin_clean),se);     
      
      % frame info
      nc_lin_indices = sub2ind(size(protein_bin_clean),nc_y_vec,nc_x_vec);
      
      % take convex hull
      stats = regionprops(nc_frame_mask,'ConvexHull');
      nc_ref_frame_raw = zeros(size(nc_frame_mask),'uint8'); 
%       zeroFilter = false(size(nc_ref_frame_raw));
      for j = 2:length(stats)
          hull_points = stats(j).ConvexHull;
          mask = poly2mask(hull_points(:,1),hull_points(:,2),yDim,xDim);   
          nc_bin_ids = mask(nc_lin_indices);
          if sum(nc_bin_ids) == 1 % enforce unique spot-nucleus-assignment
%             mask = imdilate(mask,seBig);
%             maskBig = imdilate(mask,seSmall);
            % if there are overlapping pixels, set these to zero
%             zeroFilter = zeroFilter |((logical(nc_ref_frame_raw).*mask)>0 | (maskBig&~mask));
            nc_ref_frame_raw(mask) = j;%nc_master_vec_u(nc_bin_ids);
%             nc_ref_frame_raw(zeroFilter) = 0;
%           else
%             error('wtf')
          end
      end       

%       nc_ref_frame_raw(boundarymask(nc_ref_frame_raw)) = 0;
      nc_ref_frame = nc_ref_frame_raw;%imdilate(nc_ref_frame_raw,se);
      
      nuclearMask = nc_ref_frame;

end