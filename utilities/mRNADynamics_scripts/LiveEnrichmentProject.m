classdef LiveEnrichmentProject
    %liveProject Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        Project = '';
        includedExperimentNames = [];
        includedExperiments = {};
        dataPath = '';
        figurePath = '';
        
        unhealthyNames = [];
        
        ignoredExperimentNames = [];
                
        hasSpots = [];
        hasParticles = [];
        hasSchnitzcells = [];
        hasCompiledParticles = [];
        
        anaphaseFramesAnnotated = [];
        
        hasAPInfo = [];        
        has3DSpotInfo = [];
        hasNucleusProbabilityMaps = [];
        hasProteinInfo = [];
        
        hasTracesCompiled = [];
        hasLocalProteinSampled = [];
                
    end
    
    properties (Hidden = true)
                
    end
    
    properties (Dependent = true)
        
        experimentComparisons;
        
    end
    
    methods
        %% Constructors
        
        
        function this = LiveEnrichmentProject(Project)
            %liveProject Construct an instance of this class
            %   Detailed explanation goes here
            this.Project = Project;
            [approvedPrefixes, dropboxFolder] = getProjectPrefixes(Project,'customApproved','ReadyForEnrichmentAnalysis');
            this.includedExperimentNames = string(approvedPrefixes); %NL: fiddled with this to make it specific to enrichment
            this.ignoredExperimentNames = string( getProjectPrefixes(Project,'customApproved','ReadyForEnrichmentAnalysis','inverseFlag') ); 
            
            slashes = regexp(dropboxFolder,'/|\');
            this.dataPath = [dropboxFolder(1:slashes(end)) 'ProcessedEnrichmentData' filesep Project filesep];
            this.figurePath = [dropboxFolder(1:slashes(end)) 'LocalEnrichmentFigures' filesep 'PipelineOutput' filesep Project filesep];
            
            for i = 1:length(this.includedExperimentNames)
                this.includedExperiments{i} = LiveEnrichmentExperiment(this.includedExperimentNames{i});
                isUnhealthy = this.includedExperiments{i}.isUnhealthy;
                if ~isnan(isUnhealthy) && isUnhealthy
                    this.unhealthyNames = [ this.unhealthyNames,...
                        this.includedExperimentNames{i}];
                end
            end
            
            % check progress within enrichment pipeline
            this.hasTracesCompiled = exist([this.dataPath 'spot_struct.mat'],'file');
            this.hasLocalProteinSampled = exist([this.dataPath 'spot_struct_protein.mat'],'file');
            
            
            this.hasSpots = haveSpots(this);
            this.hasParticles = haveParticles(this);
            this.hasSchnitzcells =  haveSchnitzcells(this);
            this.hasCompiledParticles = haveCompiledParticles(this);
            this.anaphaseFramesAnnotated = haveAnaphaseFrames(this);
            
            this.hasNucleusProbabilityMaps = haveNucleusProbabilityMaps(this);
            
            this.hasAPInfo = checkAPInfo(this);
            this.has3DSpotInfo = check3DSpotInfo(this);
            this.hasProteinInfo = checkProteinInfo(this);
        end
        
        
        
        
        %% Methods
        
        function hasSchnitzcells = haveSchnitzcells(this)
            
            numValidProjects = length(this.includedExperimentNames);
            
            hasSchnitzcells = zeros(1, numValidProjects);
            
            for k = 1:numValidProjects
                hasSchnitzcells(k) = ...
                    logical(this.includedExperiments{k}.hasSchnitzcellsFile);
            end
            
        end
        
        function hasParticles = haveParticles(this)
            
            numValidProjects = length(this.includedExperimentNames);
            
            hasParticles = zeros(1, numValidProjects);
            
            for k = 1:numValidProjects
                hasParticles(k) = ...
                    logical(this.includedExperiments{k}.hasParticlesFile);
            end
            
        end
        
        function hasSpots = haveSpots(this)
            
            numValidProjects = length(this.includedExperimentNames);
            
            hasSpots = zeros(1, numValidProjects);
            
            for k = 1:numValidProjects
                hasSpots(k) = ...
                    logical(this.includedExperiments{k}.hasSpotsFile);
            end
            
        end
        
        function hasCompiledParticles = haveCompiledParticles(this)
            
            numValidProjects = length(this.includedExperimentNames);
            
            hasCompiledParticles = zeros(1, numValidProjects);
            
            for k = 1:numValidProjects
                hasCompiledParticles(k) = ...
                    logical(this.includedExperiments{k}.hasCompiledParticlesFile);
            end
            
        end
        
        function hasNucleusProbabilityMaps = haveNucleusProbabilityMaps(this)
            
            numValidProjects = length(this.includedExperimentNames);
            
            hasNucleusProbabilityMaps = zeros(1, numValidProjects);
            
            for k = 1:numValidProjects
                hasNucleusProbabilityMaps(k) = ...
                    logical(this.includedExperiments{k}.hasHisProbFiles);
            end
            
        end
        
        
        function anaphaseFramesAnnotated = haveAnaphaseFrames(this)
            
            numValidProjects = length(this.includedExperimentNames);
            
            anaphaseFramesAnnotated = zeros(1, numValidProjects);
            
            for k = 1:numValidProjects
                anaphaseFramesAnnotated(k) = ...
                    nansum(this.includedExperiments{k}.anaphaseFrames) > 1;
            end
            
        end
        
        function  comparedSettings  = getExperimentComparisons(this)
            
                [comparedSettings,rawSettings] = compareExperimentSettings(this.Project);
                
        end
        
        function hasAPInfo = checkAPInfo(this)
            
            numValidProjects = length(this.includedExperimentNames);
            
            hasAPInfo = zeros(1, numValidProjects);
            
            for k = 1:numValidProjects
                APDetectionFile = [this.includedExperiments{k}.resultsFolder 'APDetection.mat'];
                hasAPInfo(k) = exist(APDetectionFile, 'file')~=0;
            end
            
        end
        
        function has3DSpotInfo = check3DSpotInfo(this)
            
            numValidProjects = length(this.includedExperimentNames);
            
            has3DSpotInfo = zeros(1, numValidProjects);
            
            for k = 1:numValidProjects
                SpotSegmetation = [this.includedExperiments{k}.resultsFolder 'Spots3DToken.mat'];
                if exist(SpotSegmetation, 'file')~=0
                    load(SpotSegmetation,'Spots3DToken');
                    has3DSpotInfo(k) = Spots3DToken >= 738110;
                end
            end
            
        end
        
        function hasProteinInfo = checkProteinInfo(this)
            
            numValidProjects = length(this.includedExperimentNames);
            
            hasProteinInfo = zeros(1, numValidProjects);
            
            for k = 1:numValidProjects                
                hasProteinInfo(k) = ~isempty(this.includedExperiments{k}.inputChannels);
            end
            
        end
        
    end
end

