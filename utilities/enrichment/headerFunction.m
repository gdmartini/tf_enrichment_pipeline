function [liveProject, numExperiments, dataName, hasAPInfo, has3DSpotInfo,...
                  hasProteinInfo, hasNucleusProbFiles] = headerFunction(projectName)

% check to see if there is a full mRNADynamics repo on the working path. If
% so remove it to prevent function cross-talk
cleanUpmRNADynamics;
       
% Determine key project characteristics
liveProject = LiveEnrichmentProject(projectName);
numExperiments = length(liveProject.includedExperimentNames);

% Make the output filepath
mkdir(liveProject.dataPath);

% Assign save names
dataName = [liveProject.dataPath 'spot_struct.mat']; % names for compiled elipse struct

% Check for optiuonal parameter felds
hasAPInfo = all(liveProject.hasAPInfo);
has3DSpotInfo = all(liveProject.has3DSpotInfo);
hasProteinInfo = all(liveProject.hasProteinInfo);
hasNucleusProbFiles = all(liveProject.hasNucleusProbabilityMaps);