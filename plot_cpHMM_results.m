% script to plot results from cpHMM inference
clear
close all
addpath(genpath('utilities'))

projectNameCell = {'EveWT'};%{'EveS1Null','EveGtSL','EveGtSL-S1Null','EveWT'};

% Set basic plotting parameters
xVar = 'fluo';

MarkerSize = 50;
blue = [115 143 193]/256;
purple = [171 133 172]/256;
red = [213 108 85]/256;

for p = 1%:length(projectNameCell)
    % set project to analyze 
    projectName = projectNameCell{p};

    % get path to results
    liveProject = LiveProject(projectName);
    resultsDir = [liveProject.dataPath filesep 'cpHMM_results' filesep];
    
    % make figure directory
    figureDir = [resultsDir 'figures' filesep];
    mkdir(figureDir);
    
    % get list of projects
    resultList = dir([resultsDir '*result*']);            
    
    for r = 1%:length(resultList)
        % load data
        load([resultsDir filesep resultList(r).name]);
      
        % get index of x axis
        dashes = strfind(resultList(r).name,'_');
        additionalVar = resultList(r).name(dashes(end)+1:strfind(resultList(r).name,'.mat')-1);
        indexVarCols = {'AP','Time','fluo', additionalVar}; %NL: this will eventually be dynamic
        xIndex = find(strcmp(indexVarCols,xVar)); 
        
        % extract index variable array
        indexVarArray = compiledResults.inferenceOptions.indexInfo.indexVarArray;
        % figure out how many unique groups we need to plot
        [newVarArray, mapTo, indexList] = unique(indexVarArray(compiledResults.groupID_index,[1:xIndex-1 xIndex+1:end]),'rows');
        % see which variables actually change
        usedGroupers = [compiledResults.inferenceOptions.apBinFlag compiledResults.inferenceOptions.timeBinFlag... 
          ~isempty(compiledResults.inferenceOptions.intensityBinVar) ~isempty(compiledResults.inferenceOptions.AdditionalGroupingVariable)];                                        
        
        lgd_flag = 0;
        if sum(usedGroupers)>2
          error('Code currently does not support plots with more than 1 additional grouping variable')
        elseif sum(usedGroupers)==2
          gpIndex = find(~ismember(1:length(indexVarCols),xIndex)&usedGroupers);
          gpVar = indexVarCols{gpIndex};
          lgd_flag = 1;
          lgd_prefix = [gpVar ' '];      
          if gpIndex > xIndex
            gpIndex = gpIndex-1;
          end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % make figures
        close all
        
        r_trend = figure;
        hm_cm = flipud(brewermap(size(newVarArray,1)+2,'Set1'));
        colormap(hm_cm);
        lgd_str = {};
        hold on
        for i = 1:size(newVarArray,1)
            ind_list = find(indexList==i);
            e = errorbar(compiledResults.fluo_mean(ind_list),compiledResults.init_vec_mean(ind_list),compiledResults.init_vec_ste(ind_list),'o','Color','black','LineWidth',1);          
            e.CapSize = 0;
        end
        
        for i = 1:size(newVarArray,1)
            ind_list = find(indexList==i);
            s(i) = scatter(compiledResults.fluo_mean(ind_list),compiledResults.init_vec_mean(ind_list),MarkerSize,'o','MarkerFaceColor',hm_cm(i+1,:),'MarkerEdgeColor','black');
            if lgd_flag
              lgd_str(i) = {[lgd_prefix num2str(newVarArray(i,gpIndex))]};
            end
        end        
        grid on
%         xlim(x_lim)
        % ylim([50 95])
        xlabel('MS2 spot intensity (au)') % NL: need to make this dynamic
        ylabel('burst amplitude (au/min)')
        if lgd_flag
            legend(s,lgd_str{:},'Location','southeast')
        end
        % set(gca,'Fontsize',14)
        StandardFigure([],gca)
        box on
%         saveas(r_trend,[FigPath,'burst_amp.tif'])
%         saveas(r_trend,[FigPath,'burst_amp.pdf'])
        %%

        dur_trend = figure;        
                
        hold on
        for i = 1:size(newVarArray,1)
            ind_list = find(indexList==i);
            e = errorbar(compiledResults.fluo_mean(ind_list),compiledResults.dur_vec_mean(ind_list),compiledResults.dur_vec_ste(ind_list),'o','Color','black','LineWidth',1);          
            e.CapSize = 0;
        end
        
        for i = 1:size(newVarArray,1)
            ind_list = find(indexList==i);
            s(i) = scatter(compiledResults.fluo_mean(ind_list),compiledResults.dur_vec_mean(ind_list),MarkerSize,'o','MarkerFaceColor',hm_cm(i+1,:),'MarkerEdgeColor','black');   
        end        
        grid on

        xlabel('MS2 spot intensity (au)') % NL: need to make this dynamic
        ylabel('burst duration (min)')
        if lgd_flag
            legend(s,lgd_str{:},'Location','southeast')
        end
        % set(gca,'Fontsize',14)
        StandardFigure([],gca)
        box on
%         saveas(dur_trend,[FigPath,'burst_dur.tif'])
%         saveas(dur_trend,[FigPath,'burst_dur.pdf'])

        %%   

        freq_trend = figure;
        
        hold on
        for i = 1:size(newVarArray,1)
            ind_list = find(indexList==i);
            e = errorbar(compiledResults.fluo_mean(ind_list),compiledResults.freq_vec_mean(ind_list),compiledResults.freq_vec_ste(ind_list),'o','Color','black','LineWidth',1);          
            e.CapSize = 0;
        end
        
        for i = 1:size(newVarArray,1)
            ind_list = find(indexList==i);
            s(i) = scatter(compiledResults.fluo_mean(ind_list),compiledResults.freq_vec_mean(ind_list),MarkerSize,'o','MarkerFaceColor',hm_cm(i+1,:),'MarkerEdgeColor','black');   
        end        
        grid on

        xlabel('MS2 spot intensity (au)') % NL: need to make this dynamic
        ylabel('burst frequency (1/min)')
        if lgd_flag
            legend(s,lgd_str{:},'Location','southeast')
        end        
        StandardFigure([],gca)
        box on
        
%         saveas(freq_trend,[FigPath,'burst_freq.tif'])
%         saveas(freq_trend,[FigPath,'burst_freq.pdf'])
    end
end