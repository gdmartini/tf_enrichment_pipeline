function spot_struct = main01_compile_traces(projectName,varargin)

%
% nucleusStruct = main01_compile_traces(dataStatusTab,dropboxFolder,varargin)
%
% DESCRIPTION
% Function to compile relevant outputs from main mRNADynamics image 
% analysis pipeline across multiple experiments
%
%
% ARGUMENTS
% projectName: master ID variable (should match a tab name in the 
%              DataStatus.xlsx spreadsheet)
%
% OPTIONS
% firstNC: script defaults to taking only nc14 traces. If you want
%          earlier traces, pass 'firstNC', followed by desired nuclear 
%          cycle number
%
% other: script allows any default variable to be set using format:
%       "VariableNameString", VariableValue
%
% OUTPUT
% spot_struct: compiled data set contain key nucleus and
%                 particle attributes, contains the following fields:
%
%     setID: ID number for this dataset (experiment), see set_key.mat to 
%            match setID to experiment Prefix
%     sourcePath: full path to the experiment's DynamicsResults folder
%     APFlag: indicates whether or not there is AP position infor for this 
%             experiment
%     threeDFlag: indicates whether or not there is 3D spot fitting data 
%                 for this experiment
%     xDim: # of pixels in X
%     yDim: # of pixels in Y
%     zDim: # of steps in Z
%     zStep: size of steps in z, in um
%     pixelSize: size of xy pixels, in um
%     particleID: 
%     xPosParticle: 
%     yPosParticle: 
%     zPosParticle: 
%     APPosParticle: 
%     spotFrames: 
%     xPosParticle3D: 
%     yPosParticle3D: 
%     zPosParticle3D: 
%     fluo3D: 
%     fluo: 
%     fluoOffset: 
%     qcFlag: 
%     N: ???
%     sparsity: 
%     xPosNucleus: x position of center of nucleus
%     yPosNucleus: y position of center of nucleus
%     rawNCPprotein: 
%     frames: 
%     nucleusID
%     ncID
%     ncStart
%     minDP
%     minTime
%     pctSparsity
% 
%     time
% 
%     sisterIndex
%     sisterParticleID

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% Set defaults %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all force
addpath(genpath('utilities'));
% Defaults
firstNC = 14;   % first nuclear cycle to pull data from
minDP = 14;     % what is this for?
pctSparsity = 50;   % 
twoSpotFlag = contains(projectName, '2spot');
minTime = 0*60; % take no fluorescence data prior to this point
tresInterp = 20; 
calculatePSF = false;
% projectName = projectName;
SpotChannelIndex = 1; % this does nothing at the moment, but can serve as a starting point if ever we wish to analyze two-spot-two-color data 
ncRefVec = 9:14;

%% %%%%%%%%%%%%%%%%%%%%%%% Process input parameters %%%%%%%%%%%%%%%%%%%%%%%
includeFlaggingInfo = false; % Added by GM 12/3/20
for i = 1:numel(varargin)
    if ischar(varargin{i})
        if strcmpi(varargin{i}, 'firstNC')
            if  i < numel(varargin) && isnumeric(varargin{i+1}) && (9 <= varargin{i+1} && varargin{i+1} <= 14)
                firstNC = varargin{i+1};
            else
                error(['The ''firstNC'' option must be followed by an integer from 9 to 14, inclusive. Default firstNC is ' num2str(firstNC) '.']) 
            end
        elseif strcmpi(varargin{i}, 'projectName')
            error('Pipeline currently does not support alternate names for projects (otehr than what is on the data status tab)')
        elseif strcmpi(varargin{i}, 'includeFlaggingInfo')
            includeFlaggingInfo = true; % GM 12/3/20: Functionality for this option has yet to be implemented
        elseif i < numel(varargin)                         
            eval([varargin{i} '=varargin{i+1};']);        
        end
    end
end

%% %%%%%%%%%%%%%%%%%%%%% Set Path Specs, ID Vars %%%%%%%%%%%%%%%%%%%%%%%%%%

[liveProject, numExperiments, dataName, hasAPInfo, has3DSpotInfo, hasProteinInfo] = headerFunction(projectName);

%% %%%%%%%%%%%%%%%%% Extract relevant processed data %%%%%%%%%%%%%%%%%%%%%%

% Generate master structure with info on all nuclei and traces in
% constituent sets
spot_struct = [];

% Add data from each experiment to the master structure
h = waitbar(0,'Compiling data ...');
for i = 1:numExperiments
    
    waitbar(i/numExperiments,h, ['Compiling data for dataset ' num2str(i) ' of ' num2str(numExperiments)])
    
    setID = i;
    currExperiment = liveProject.includedExperiments{i};
        
    %%%%%%%% Read in processed data from main mRNADyanmics pipeline %%%%%%%
    
    schnitzcells = getSchnitzcells(currExperiment);
    processedData = getCompiledParticles(currExperiment);   % this structure contains all info compiled for this experiment
    compiledParticles = processedData.CompiledParticles;    % this is the cell array containing the CompiledParticles cell array
    if iscell(compiledParticles)
        compiledParticles = compiledParticles{SpotChannelIndex};               % this assumes there is only one channel with spots to analyze
    end   
    
    % MT 2020-08-11: Right now we don't use any of the data in
    % CompiledNuclei (from mRNADynamics pipeline) because we generate all
    % of it ourselves in main02. This might change in the future, so I'm
    % including this section so it's easy to use should we need it     
%     compiledNuclei = processedData.CompiledNuclei;        % this is the CompiledNuclei structure we actually want to use  
%     if iscell(compiledNuclei)
%         compiledNuclei = compiledNuclei{1};                     % unlike compiledParticles, I don't think this gets stored in a cell array, but just in case...           
%     end
    
    
    %%%%%%%%%%%%%%%%%%% Get basic data characteristics %%%%%%%%%%%%%%%%%%%%
        
    % Check if certain info is present for this experiment        
    DetrendedZFlag = isfield(compiledParticles, 'zPosDetrended');   %Z pos that accounts for stack repositioning

    % Pull trace, time, and frame variables
    timeRaw = processedData.ElapsedTime*60; %time vector, [sec] 
    framesRaw = 1:length(timeRaw); % Frame list
    tracesRaw = processedData.AllTracesVector;  %array with a column for each trace
    if iscell(tracesRaw)        %check to see if traces are stored in cell array
        tracesRaw = tracesRaw{SpotChannelIndex};   
    end
    
    % Find start frames for nc's we wish to include. Pull this from
    % FrameInfo if not in movie database
    anaphaseFrames = currExperiment.anaphaseFrames;
    if all(isnan(anaphaseFrames))
      FrameInfo = getFrameInfo(currExperiment);
      ncVec = [FrameInfo.nc];
      ncIndex = unique(ncVec);
      for n = 1:length(ncIndex)
        anaphaseFrames(ncRefVec==ncIndex(n)) = find(ncVec==ncIndex(n),1);
      end
    end
    ncIndices = find(ncRefVec>=firstNC&~isnan(anaphaseFrames'));
    ncStartFrameVec = anaphaseFrames(ncIndices)';
    ncStartFrameVec(end) = max([1 ncStartFrameVec(end)]); 
    ncStartFrameVec(end+1) = framesRaw(end)+1;
    firstTime = timeRaw(min(ncStartFrameVec));            
    
    % iterate through nuclear cycles
    for ncInd = 1:length(ncStartFrameVec)-1
      
        % initialize structure to store nucleus and particle info
        compiledSchnitzCells = struct;

        nucleusCounter = 1;    %counter to ensure no empty spaces in the compiledSchnitzCells
                        %struct if some schnitzcells have no frames in the
                        %desired nuclear cycle
      
        firstNCFrame = ncStartFrameVec(ncInd);%max([1, processedData.(['nc' num2str(firstNC)])]);      
        lastNCFrame = ncStartFrameVec(ncInd+1)-1;

        % filter reference vectors
        framesNC = framesRaw(firstNCFrame:lastNCFrame); 
        tracesNC = tracesRaw(firstNCFrame:lastNCFrame,:);
        timeNC = timeRaw(firstNCFrame:lastNCFrame);    
        timeNC = timeNC - firstTime; %normalize to start of first nc

        %%%%%%%%%%%%%%%%%%%% Compile nucleus schnitz info %%%%%%%%%%%%%%%%%%%%%
      
        for s = 1:length(schnitzcells)
            schnitzFrames = schnitzcells(s).frames;

            % Filter for only the frames in this nc
            ncFilter = ismember(schnitzFrames,framesNC);
            rawNucleusFrames = schnitzFrames(ncFilter); 

            if length(rawNucleusFrames) >= 1     %only grab data from the desired nuclear cycle(s)

                % Add info that's the same for all schnitzcells in this
                % experiment
                compiledSchnitzCells(nucleusCounter).setID = setID;                                                
                compiledSchnitzCells(nucleusCounter).nc = ncRefVec(ncIndices(ncInd));

                % Initialize particle fields--will be set to particle values 
                % for nuclei with matching particle
                compiledSchnitzCells = initializeParticleFields(...
                  compiledSchnitzCells,nucleusCounter,sum(ncFilter),has3DSpotInfo,hasAPInfo);              
                
                % Add QC-related flags
                compiledSchnitzCells(nucleusCounter).qcFlag = NaN;
                compiledSchnitzCells(nucleusCounter).N = NaN;
                compiledSchnitzCells(nucleusCounter).sparsity = NaN;

                % Add core nucleus info
                x = schnitzcells(s).cenx;            
                y = schnitzcells(s).ceny;
                compiledSchnitzCells(nucleusCounter).xPosNucleus = x(ncFilter);
                compiledSchnitzCells(nucleusCounter).yPosNucleus = y(ncFilter);

                % Add protein and nucleus info
                if hasProteinInfo
                  compiledSchnitzCells(nucleusCounter).rawNCPprotein = nanmax(schnitzcells(s).Fluo(ncFilter,:),[],2);
                end
                compiledSchnitzCells(nucleusCounter).frames = rawNucleusFrames';            
                compiledSchnitzCells(nucleusCounter).nucleusID = s;     
                compiledSchnitzCells(nucleusCounter).ncID = eval([num2str(setID) '.' sprintf('%04d',s)]);  
                compiledSchnitzCells(nucleusCounter).minDP = minDP;
                compiledSchnitzCells(nucleusCounter).minTime = minTime;
                compiledSchnitzCells(nucleusCounter).pctSparsity = pctSparsity;

                % Add time and set info
                compiledSchnitzCells(nucleusCounter).time = timeNC(ismember(framesNC,rawNucleusFrames));
                if numel(unique(compiledSchnitzCells(nucleusCounter).time))~=numel(compiledSchnitzCells(nucleusCounter).time)
                    error(['Non-unique time values. Check FrameInfo for Prefix: ' currExperiment.Prefix])
                end
                % initialize nucleus qc flag
                compiledSchnitzCells(nucleusCounter).missingNucleusFrames = 0; 
                
                % Add sister spot fields
                compiledSchnitzCells(nucleusCounter).sisterIndex = NaN;
                compiledSchnitzCells(nucleusCounter).sisterParticleID = NaN; 
            
                % increment counter to ensure no empty spaces in the 
                % compiledSchnitzCells struct if some schnitzcells have no 
                % frames in the desired nuclear cycle
                nucleusCounter = nucleusCounter + 1;               
            end
        end
    
        % Index vector to cross-ref w/ particles            
        schnitzIndex = [compiledSchnitzCells.nucleusID];          
        
        % Iterate through traces
        for j = 1:size(tracesNC,2)  
            % Raw fluo trace
            rawTrace = tracesNC(:,j); 

            % Get nucleus ID
            schnitz = compiledParticles(j).schnitz;    

            % Find first and last expression frames, requiring that spots
            % cannot apear earlier than some fixed time (minTime)
            traceStart = find(timeNC>=minTime&~isnan(rawTrace'),1);
            traceStop = find(~isnan(rawTrace),1,'last');        

            % skip particles not in specified nc range
            if isempty(traceStart)
                continue
            end   

            %Create versions with all intervening frames present (missing frames
            %appear as NaNs)
            traceFull = rawTrace(traceStart:traceStop)';                                   
            traceFramesFull = framesNC(traceStart:traceStop);          

            % Perform qc tests    
            nDP = sum(~isnan(traceFull));
            sparsity = prctile(diff(find(~isnan(traceFull))),pctSparsity);
            qcFlag = nDP >= minDP && sparsity == 1;     

            % trace-nucleus mapping may be many-to-1
            ncIndex = find(schnitzIndex==schnitz);  
            if length(ncIndex) ~= 1
                error(['Problem with particle-nucleus crossreference for Prefix: ' currExperiment.Prefix])          
            end 

            % Identifier variable                        
            particle = compiledParticles(j).OriginalParticle;                        
            particleID = eval([num2str(setID) '.' sprintf('%04d',particle)]);

            % check to see if a different particle has already been assigned
            % if so, create a new entry
            if ~isnan(compiledSchnitzCells(ncIndex).particleID) && twoSpotFlag            
                compiledSchnitzCells(end+1) = compiledSchnitzCells(ncIndex);

                % update cross-reference variables
                compiledSchnitzCells(ncIndex).sisterParticleID = particleID;
                compiledSchnitzCells(end).sisterParticleID = compiledSchnitzCells(ncIndex).sisterParticleID;
                compiledSchnitzCells(ncIndex).sisterIndex = length(compiledSchnitzCells);
                compiledSchnitzCells(end).sisterIndex = ncIndex;

                % redefine index variables
                ncIndex = length(compiledSchnitzCells);

                % reset particle fields
                nEntries = length(compiledSchnitzCells(ncIndex).xPosParticle);            
                compiledSchnitzCells(ncIndex).ParticleID = NaN;

                % Initialize particle fields--will be set to particle values 
                % for nuclei with matching particle
                compiledSchnitzCells = initializeParticleFields(...
                  compiledSchnitzCells,ncIndex,nEntries);%,threeDFlag,APFlag,DetrendedZFlag);           
            end

            % Record particle identifier
            compiledSchnitzCells(ncIndex).particleID = particleID;

            % Find overlap between nucleus and trace
            rawNucleusFrames = compiledSchnitzCells(ncIndex).frames;         
            spotFilter = ismember(rawNucleusFrames,traceFramesFull); 
            
            compiledSchnitzCells(ncIndex).spotFrames = spotFilter;            
            if sum(spotFilter) < length(traceFramesFull)            
                warning(['Inconsistent particle and nucleus frames for Prefix: ' currExperiment.Prefix ' .Flagging...'])  
                compiledSchnitzCells(ncIndex).missingNucleusFrames = 1;
            end
            
            % record fluorescence info             
            compiledSchnitzCells(ncIndex).fluo(spotFilter) = traceFull(ismember(traceFramesFull,rawNucleusFrames)); 
                        
            % Find intersection btw full frame range and CP frames        
            rawParticleFrames = compiledParticles(j).Frame;
%             rawParticleFrames = rawParticleFrames(ismember(rawParticleFrames,traceFramesFull));
            
            % make filters
            ncSpotFilter1 = ismember(rawNucleusFrames,rawParticleFrames);
            ncSpotFilter2 = ismember(rawParticleFrames,rawNucleusFrames);
                    
            % add offset info
            compiledSchnitzCells(ncIndex).fluoOffset(ncSpotFilter1) = compiledParticles(j).Off(ncSpotFilter2);
            % x, y, and z info                                
            compiledSchnitzCells(ncIndex).xPosParticle(ncSpotFilter1) = compiledParticles(j).xPos(ncSpotFilter2);
            compiledSchnitzCells(ncIndex).yPosParticle(ncSpotFilter1) = compiledParticles(j).yPos(ncSpotFilter2);   
            compiledSchnitzCells(ncIndex).zPosParticle(ncSpotFilter1) = compiledParticles(j).zPos(ncSpotFilter2);
            if hasAPInfo
                compiledSchnitzCells(ncIndex).APPosParticle(ncSpotFilter1) = compiledParticles(j).APposParticle(ncSpotFilter2)*100;
            end
            % 3D info
            if has3DSpotInfo
                compiledSchnitzCells(ncIndex).xPosParticle3D(ncSpotFilter1) = compiledParticles(j).xPosGauss3D(ncSpotFilter2);            
                compiledSchnitzCells(ncIndex).yPosParticle3D(ncSpotFilter1) = compiledParticles(j).yPosGauss3D(ncSpotFilter2);              
                compiledSchnitzCells(ncIndex).zPosParticle3D(ncSpotFilter1) = compiledParticles(j).zPosGauss3D(ncSpotFilter2);            
                compiledSchnitzCells(ncIndex).fluo3D(ncSpotFilter1) = compiledParticles(j).Fluo3DRaw(ncSpotFilter2);
            end
            if DetrendedZFlag
                compiledSchnitzCells(ncIndex).zPosParticleDetrended(ncSpotFilter1) = compiledParticles(j).zPosDetrended(ncSpotFilter2);            
            end            
            % add qc info
            compiledSchnitzCells(ncIndex).N = nDP;
            compiledSchnitzCells(ncIndex).sparsity = sparsity;        
            compiledSchnitzCells(ncIndex).qcFlag = qcFlag;  
        end      
        spot_struct = [spot_struct  compiledSchnitzCells];
    end
end
close(h)

% add additional fields
% NL: this does nothing at the moment, but leaving to preserve two-spot
% compatibility
for i = 1:numel(spot_struct)
    spot_struct(i).twoSpotFlag = twoSpotFlag;
    spot_struct(i).targetLocusFlag = NaN;
    spot_struct(i).controlLocusFlag = NaN;
end    

disp('interpolating data...')

% generate interpolation fields
interpFields = {'fluo'};
if has3DSpotInfo
    interpFields(end+1) = {'fluo3D'};
end
if hasAPInfo
    interpFields(end+1) = {'APPosParticle'};
end
if hasProteinInfo
    interpFields(end+1) = {'rawNCPprotein'};
end

% calculate interpolation time
interpGrid = 0:tresInterp:60*60;

% use 95th percentile of fluorescence values to set scale for flagging 
% unusual jump events
fluo_scale = prctile([spot_struct.fluo],97.5);
big_blip_thresh = fluo_scale;

for i = 1:numel(spot_struct)
    timeVec = spot_struct(i).time;
    fluoVec = spot_struct(i).fluo;
    startIndex = find(~isnan(fluoVec),1);
    stopIndex = find(~isnan(fluoVec),1,'last');    
    timeVec = timeVec(startIndex:stopIndex);       
    if length(timeVec)>1
        timeInterp = interpGrid(interpGrid>=timeVec(1)&interpGrid<=timeVec(end));
        for  j = 1:numel(interpFields)
            vec = spot_struct(i).(interpFields{j})(startIndex:stopIndex);
            
            if contains(interpFields{j},'fluo')
                %Look for clusters of 6 or more NaNs
                kernel = [1,1,1,1,1];
                vecNans = isnan(vec);
                vecConv = conv(kernel,vecNans);
                vecConv = vecConv(3:end-2);
                zIDs = find(vecConv==numel(kernel));
                zIDs = unique([zIDs-1 zIDs zIDs+1]); % get set of z_ids    
                vec(zIDs) = 0; % set clusters to zeros    
                vec(vec<0) = 0; % deal with negative values    

                % find single dp "blips". These will be replaced via interpolation
                % interpolate remaining NaNs  
                fluo_dd = [0 diff(vec,2) 0];
                vec(abs(fluo_dd)>big_blip_thresh) = NaN;

                % find and remove suspsiciously large rises
                fluo_d = [0 diff(vec,1) 0];
                vec(abs(fluo_d)>big_blip_thresh*0.75) = NaN;
            end
            
            % interpolate remaining NaNs    
            queryPoints = timeVec(isnan(vec));
            referenceTime = timeVec(~isnan(vec));
            referenceVec = vec(~isnan(vec));
%             if ~isempty(queryPoints)
%                 newF = interp1(referenceTime,referenceVec,queryPoints);  % NL: this is a little weird but leaving for now
%                 vec(ismember(timeVec,queryPoints)) = newF;   
%             else
%                 vec = referenceVec;
%             end                 
            if length(referenceTime)>1
                % Interpolate to standardize spacing        
                spot_struct(i).([interpFields{j} 'Interp']) = interp1(referenceTime,referenceVec,timeInterp);
            elseif length(referenceTime)==1
                spot_struct(i).([interpFields{j} 'Interp']) = referenceVec;
            else
                spot_struct(i).([interpFields{j} 'Interp']) = NaN(size(timeInterp));
            end              
        end
    else
        timeInterp = NaN;
        for  j = 1:numel(interpFields)
            spot_struct(i).([interpFields{j} 'Interp']) = NaN;
        end
    end
    spot_struct(i).timeInterp = timeInterp;
    spot_struct(i).tresInterp = tresInterp;
end

% save
save(dataName ,'spot_struct') 

if calculatePSF % NL: this needs to be fixed...
    disp('calculating psf dims...')
    % call function to calculate average psf difs
    calculate_average_psf(projectName,dropboxFolder);
end

disp('done.')
