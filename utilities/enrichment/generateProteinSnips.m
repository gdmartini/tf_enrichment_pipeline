function [spot_protein_snips_mixed, spot_mcp_snips_mixed, edge_control_protein_snips_mixed, edge_control_mcp_snips_mixed] = ...
                                      generateProteinSnips(snip_data, DistFilter)

% extract snip stacks
spot_protein_snips = cat(3,snip_data.spot_protein_snips);
spot_mcp_snips = cat(3,snip_data.spot_mcp_snips);
edge_control_protein_snips = cat(3,snip_data.edge_control_protein_snips);
edge_control_mcp_snips = cat(3,snip_data.edge_control_mcp_snips);

% apply distance filter
spot_protein_snips_dist = spot_protein_snips(:,:,DistFilter);
edge_control_protein_snips_dist = edge_control_protein_snips(:,:,DistFilter);
spot_mcp_snips_dist = spot_mcp_snips(:,:,DistFilter);
edge_control_mcp_snips_dist = edge_control_mcp_snips(:,:,DistFilter);

% randomize snip orientation
snip_size = size(spot_mcp_snips_dist,1);
inv_mat = [fliplr(1:snip_size); 1:snip_size]' ;
spot_protein_snips_mixed = NaN(size(spot_protein_snips_dist));
edge_control_protein_snips_mixed = NaN(size(spot_protein_snips_dist));
spot_mcp_snips_mixed = NaN(size(spot_protein_snips_dist));
edge_control_mcp_snips_mixed = NaN(size(spot_protein_snips_dist));
for i = 1:size(spot_protein_snips_dist,3)
    h = 1;%ceil(rand()*2);
    v = 1;%ceil(rand()*2);
    
    spot_protein_snips_mixed(:,:,i) = spot_protein_snips_dist(inv_mat(:,v),inv_mat(:,h),i);
    edge_control_protein_snips_mixed(:,:,i) = edge_control_protein_snips_dist(inv_mat(:,v),inv_mat(:,h),i);
    spot_mcp_snips_mixed(:,:,i) = spot_mcp_snips_dist(inv_mat(:,v),inv_mat(:,h),i);
    edge_control_mcp_snips_mixed(:,:,i) = edge_control_mcp_snips_dist(inv_mat(:,v),inv_mat(:,h),i);
end