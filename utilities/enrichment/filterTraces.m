function [trace_struct_filtered, indexInfo, inferenceOptions] = filterTraces(inferenceOptions,analysis_traces)
  % this function filters traces and generates grouping variables as
  % appropriate
  % Supported grouping variables include: time, APPositionParticle, and nuclear_protein_vec
    
  % apply QC filter calculated in main01
  analysis_traces = analysis_traces([analysis_traces.qcFlag]==1);
  
  if inferenceOptions.ProteinBinFlag && inferenceOptions.FluoBinFlag  
    error('Protein and spot fluorescence binning are mutually exclusive inference options')
  end
  
  % generate stripped down data structure
  trace_struct_filtered = struct;
  
  i_pass = 1;
  for i = 1:length(analysis_traces)    
      if inferenceOptions.fluo3DFlag
          fluo = analysis_traces(i).fluo3DInterp;
      else
          fluo = analysis_traces(i).fluoInterp;
      end
      time = analysis_traces(i).timeInterp;    
      if ~isfield(inferenceOptions,'Tres')
        inferenceOptions.Tres = time(2)-time(1);
      end
      time_raw = analysis_traces(i).time;
      
      if isfield(analysis_traces,'APPosParticle')
        ap_raw = analysis_traces(i).APPosParticle;
      else
        ap_raw = ones(size(time));
      end
      
      for a = 1:length(inferenceOptions.apBins)-1 % Note: if either option is no activated, there will be only one bin
          apBounds = inferenceOptions.apBins(a:a+1);

          for t = 1:length(inferenceOptions.timeBins)-1
            
              start_time = time_raw(1);
              if inferenceOptions.truncInference(t)
                start_time = start_time + (1+inferenceOptions.nSteps)*inferenceOptions.Tres;
              end

              timeBounds = inferenceOptions.timeBins(t:t+1);

              ap_time_filter_raw = time_raw >= timeBounds(1) & time_raw < timeBounds(2) & time_raw >= start_time &...
                               ap_raw >= apBounds(1) & ap_raw < apBounds(2);

              nRaw = sum(ap_time_filter_raw);

              if nRaw >= inferenceOptions.minDP
                ap_interp = interp1(time_raw,ap_raw,time);
%                 time_vec_temp = time_raw(ap_time_filter_raw);
                
                ap_time_filter_interp = time >= timeBounds(1) & time < timeBounds(2) & time >= start_time &...
                               ap_interp >= apBounds(1) & ap_interp < apBounds(2);
%                 ap_time_filter_interp = time>=time_vec_temp(1) & time<=time_vec_temp(end);
         
                % generate new entry    
                trace_struct_filtered(i_pass).fluo = fluo(ap_time_filter_interp);
                trace_struct_filtered(i_pass).time = time(ap_time_filter_interp);
                if inferenceOptions.ProteinBinFlag || inferenceOptions.FluoBinFlag          
                    vec = analysis_traces(i).(inferenceOptions.intensityBinVar);
                    trace_struct_filtered(i_pass).mean_intensity = nanmean(vec(ap_time_filter_raw));
                end       
                trace_struct_filtered(i_pass).particleID = analysis_traces(i).particleID;  
                trace_struct_filtered(i_pass).N = sum(ap_time_filter_interp);    
                trace_struct_filtered(i_pass).apBin = a;
                trace_struct_filtered(i_pass).timeBin = t;
                trace_struct_filtered(i_pass).interp_filter = ap_time_filter_interp;
                trace_struct_filtered(i_pass).raw_filter = ap_time_filter_raw;
                if ~isempty(inferenceOptions.AdditionalGroupingVariable)
                  traceVarName = inferenceOptions.AdditionalGroupingVariable;
                  trace_struct_filtered(i_pass).(traceVarName) = analysis_traces(i).(traceVarName);
                end
                % increment
                i_pass = i_pass + 1;
              end
          end
      end
  end
  
  % make new indexing vectors
  ap_group_vec = [trace_struct_filtered.apBin];
  time_group_vec = [trace_struct_filtered.timeBin];
  % initialize dummy grouping variable  
  additional_group_vec = ones(size(ap_group_vec));
  
  % check for additional binning variable (only 1 addtional variable
  % supported for now)  
  if ~isempty(inferenceOptions.AdditionalGroupingVariable)
    % find corresponding grouper variable
    traceVarName = inferenceOptions.AdditionalGroupingVariable;
    additional_group_vec = [trace_struct_filtered.(traceVarName)];
    inferenceOptions.additionalGroupIDs = unique(additional_group_vec(~isnan(additional_group_vec)));
  end
  
  nan_filter1 = isnan(ap_group_vec) | isnan(time_group_vec) | isnan(additional_group_vec);
  % perform protein binning if appropriate
  if inferenceOptions.ProteinBinFlag || inferenceOptions.FluoBinFlag  
    
    for a = 1:length(inferenceOptions.apBins)-1 % Note: if either option is no activated, there will be only one bin                
        for t = 1:length(inferenceOptions.timeBins)-1          
            for g = 1:length(inferenceOptions.additionalGroupIDs)
                group_ids = find(ap_group_vec==a & time_group_vec==t & additional_group_vec == inferenceOptions.additionalGroupIDs(g));          
                % estimate number of bins 
                if inferenceOptions.automaticBinning
                    nTotal = 0.98*sum([trace_struct_filtered(group_ids).N]);
                    n_intensity_bins = ceil(nTotal/inferenceOptions.SampleSize);
                end

                % generate list of average protein levels
                intensity_list = [trace_struct_filtered(group_ids).mean_intensity];
                % generate protein groupings    
                q_vec = linspace(.1,.99,n_intensity_bins+1);        
                intensity_prctile_vec = quantile(intensity_list,q_vec);    
                
                % assign traces to groups    
                id_vec = discretize(intensity_list,intensity_prctile_vec);

                for i = 1:length(id_vec)
                    trace_struct_filtered(group_ids(i)).intensity_bin = id_vec(i);
                    trace_struct_filtered(group_ids(i)).intensityquantiles = intensity_prctile_vec;
                end
            end
        end
        % fill in traces that do not fit in any "Additional group" with
        % NaNs
        for i = find(nan_filter1)
          trace_struct_filtered(i).intensity_bin = NaN;
          trace_struct_filtered(i).intensity_quantiles = [-NaN NaN];
        end
    end
  else
      for i = 1:length(trace_struct_filtered)
          trace_struct_filtered(i).intensity_bin = 1;
          trace_struct_filtered(i).intensity_quantiles = [-Inf Inf];
      end
  end
  
  
  % remove traces where one more more ID field in NAN
  intensity_group_vec = [trace_struct_filtered.intensity_bin];
  
  nan_filter = isnan(intensity_group_vec) | nan_filter1;
  trace_struct_filtered = trace_struct_filtered(~nan_filter);
  ap_group_vec = ap_group_vec(~nan_filter);
  time_group_vec = time_group_vec(~nan_filter);
  intensity_group_vec = intensity_group_vec(~nan_filter);
  additional_group_vec = additional_group_vec(~nan_filter);
  
  % generate indexing structure
  indexInfo = struct;  
  [indexInfo.indexVarArray, mapTo, indexInfo.indexList] = unique([ap_group_vec' time_group_vec' intensity_group_vec' additional_group_vec'],'rows');
  indexInfo.indexVecUnique = 1:size(indexInfo.indexVarArray,1);
  indexInfo.ap_group_vec = ap_group_vec(mapTo);
  indexInfo.time_group_vec = time_group_vec(mapTo);  
  indexInfo.intensity_group_vec = intensity_group_vec(mapTo);
  indexInfo.additional_group_vec = additional_group_vec(mapTo);
  
  % update sample size
  inferenceOptions.SampleSize = repelem(inferenceOptions.SampleSize,size(indexInfo.indexVarArray,1));
  