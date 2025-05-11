function neighbor_list = find_neighbor_voxel(voxel_i, dist_voxel, pos_voxel)
% find neighbor in cortex
% Author: Jinlong Li and Guoyuan Yang, BIT.
% 
% Inputs:
%     - voxel_i
%       index of voxel in fs_LR_32k
% 
%     - dist_voxel
%       threshold of distance between boxel
%     - pos_voxel
%       no use
% 
% Outputs:
%     - neighbor_list
%       result of neighbor list    

    neighbor_list = [];
    voxel_d = 2;
    % dist < 2mm
    idx = find(dist_voxel(voxel_i, :) <= voxel_d);
    neighbor_list = [neighbor_list idx];

    % share a corner
%     pos_varia_temp = [1 1 1; 1 1 -1; 1 -1 1; 1 -1 -1; -1 1 1; -1 1 -1; -1 -1 1; -1 -1 -1];
%     pos_varia_temp = [1 1 1; -1 -1 -1];
%     pos_temp = pos_voxel(voxel_i) + pos_varia_temp .* voxel_d ;
%     for cornor_i = 1:size(pos_temp, 1)
%         [rowIdx, ~] = ismember(pos_voxel, pos_temp(cornor_i, :), 'rows');
%         neighbor_list = [neighbor_list find(rowIdx ~= 0)];
%     end
    neighbor_list = unique(neighbor_list);
end