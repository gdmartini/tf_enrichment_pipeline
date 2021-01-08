function [spot_struct_protein, new_vec_fields] = initializeProteinFields(spot_struct_protein, has3DSpotInfo)

new_vec_fields = {'spot_protein_vec', ... % protein levels at locus                 
                  'serial_null_protein_vec',... % protein levels at simulated spot                              
                  'edge_null_protein_vec',... % protein level at spot selected to be equidistant from nucleus edge                  
                  'nuclear_protein_vec',... % average inside nucleus
                  'spot_mcp_vec',... % MCP channel intensity at spot
                  'serial_null_mcp_vec',...
                  'edge_null_mcp_vec',... % " "
                  'serial_qc_flag_vec',... % Flags when there are qc concerns for virtual spot
                  'edge_qc_flag_vec', ... % " "
                  'edge_null_x_vec', ... % x position for edge-controled simulated spot
                  'edge_null_y_vec', ... " "
                  'serial_null_x_vec',... " "
                  'serial_null_y_vec',... " "              
                  'edge_null_nc_vec',... % Nucleus where control spot was sampled from (not always same as spot)
                  'spot_edge_dist_vec',... % Spot distance from nuclear boundary
                  'serial_null_edge_dist_vec'}; % Virtual spot distance from edge
% if has3DSpotInfo
%   additionalFields =   {'spot_protein_vec_3d',...                                    
%                         'serial_null_protein_vec_3d',...                  
%                         'edge_null_protein_vec_3d',... 
%                         'nuclear_protein_vec_3d',...
%                         };
%   new_vec_fields = [new_vec_fields additionalFields];
% end

% Initialize fields
for i = 1:length(spot_struct_protein)
    refVec = spot_struct_protein(i).xPosParticle;
    for j = 1:length(new_vec_fields)
        spot_struct_protein(i).(new_vec_fields{j}) = NaN(size(refVec));
    end
end