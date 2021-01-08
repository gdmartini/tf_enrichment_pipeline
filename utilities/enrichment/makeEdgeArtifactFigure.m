function makeEdgeArtifactFigure()

% First look for presence of edge artifact
dist_sigma = 0.1; %(um)
protein_delta_vec =  spot_protein_vec - null_protein_vec;
dist_index = 0:0.1:floor(prctile(dist_vec,99)*10)/10;

delta_dist_mat = NaN(numel(dist_index),NBoots);
null_dist_mat = NaN(numel(dist_index),NBoots);
spot_dist_mat = NaN(numel(dist_index),NBoots);

for n = 1:NBoots
    s_ids = randsample(1:numel(protein_delta_vec),numel(protein_delta_vec),true);
    dv_samp1 = protein_delta_vec(s_ids);    
    nn_samp = spot_protein_vec(s_ids);
    sp_samp = null_protein_vec(s_ids);
    dist_samp = dist_vec(s_ids);        
    for t = 1:numel(dist_index)
        d_weights = exp(-.5*((dist_samp-dist_index(t))/dist_sigma).^2);
        delta_dist_mat(t,n) = nansum(dv_samp1.*d_weights) / nansum(d_weights);
        null_dist_mat(t,n) = (nansum(nn_samp.*d_weights) / nansum(d_weights));
        spot_dist_mat(t,n) = (nansum(sp_samp.*d_weights) / nansum(d_weights));
    end
end
% calcualte mean and standard error
delta_dist_mean = nanmean(delta_dist_mat,2);
delta_dist_ste = nanstd(delta_dist_mat,[],2);

null_dist_mean = nanmean(null_dist_mat,2);

%% make dist-dependent fold enrichment figure and select distance threshold
pass = 0;
while ~pass
    delta_dist_fig = figure;
    e = errorbar(dist_index,100* delta_dist_mean ./ null_dist_mean,100 * delta_dist_ste ./ null_dist_mean);
    e.CapSize = 0;
    hold on
    xlabel('distance from edge (\mu m)')
    ylabel('apparent % enrichment')
    
    grid on
    y_max = 100 * nanmax(delta_dist_mean ./ null_dist_mean);
    y_min = 100 * nanmin(delta_dist_mean ./ null_dist_mean);
    p = plot([DistLim DistLim],[y_min y_max],'Color', 'red');
    legend(p, 'current limit')
    ylim([min([0 , y_min-.1*y_min*sign(y_min)]) 1.1*y_max])  
    if ManualDistThreshold
        title(strvcat('Setting Edge Distance Threshold:',...
                'If current value is staisfactory press "Enter"',...
                'Else click a desired cutoff and click "Enter"'))
        [x,~] = ginput;
        if isempty(x)
            pass = 1;
        end
    else
        pass = 1;
    end
end
title(['Enrichment vs. Distance from Nucleus Edge (' id_string ')'])
saveas(delta_dist_fig, [FigurePath write_string '_edge_artifact_plot.png'])