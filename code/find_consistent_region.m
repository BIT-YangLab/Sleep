function label_mp = find_consistent_region(labels, mesh)
% YLAB_find_local_consistent_region(labels, mesh) to find all consistent region in a brain atlas
% labels is 
%   mesh is 
%
%
%   Examples
%     YLAB_find_local_consistent_region(label, 'fs_LR_32k');
%
%   Copyright 
% find all consistent region in a brain atlas
% Author: Xinyu Wu and Guoyuan Yang, BIT.
% 
% Inputs:
%     - labels
%       the brain parcellation with each vertice or voxel matching one network
% 
%     - mesh
%       the surface template, like fs_LR_32k, fsaverageX
%
% Outputs:
%     - label_mp
%       result label


if contains(mesh, 'fs_LR')
    lh_targ_mesh = CBIG_read_fslr_surface('lh', mesh,'inflated','medialwall.annot');
    rh_targ_mesh = CBIG_read_fslr_surface('rh', mesh,'inflated','medialwall.annot');
    no_medial_wall_index = [find(lh_targ_mesh.MARS_label == 2); find(rh_targ_mesh.MARS_label == 2) + length(lh_targ_mesh.MARS_label)];
    lh_end_idx = 32492;
elseif contains(mesh, 'fsaverage') 
    lh_targ_mesh = CBIG_ReadNCAvgMesh('lh', mesh, 'inflated', 'cortex');
    rh_targ_mesh = CBIG_ReadNCAvgMesh('rh', mesh, 'inflated', 'cortex');
    
    no_medial_wall_index = [find(lh_targ_mesh.MARS_label == 2) find(rh_targ_mesh.MARS_label == 2) + length(lh_targ_mesh.MARS_label)]';
    lh_end_idx = 40962;
end


    lh_vn = lh_targ_mesh.vertexNbors ;
    rh_vn = rh_targ_mesh.vertexNbors ;
    rh_vn(rh_vn ~= 0) = rh_vn(rh_vn ~= 0) + lh_end_idx;
    neighbors = [lh_vn rh_vn];

% some area may be non, seen as one region

label_mp = zeros(64984,1);
isvisit_mp = ones(64984, 1);
isvisit_mp(no_medial_wall_index) = 0;


label_cnt = 1;
ver_list = [];
    
% bfs
for j = 1:length(labels)
    if isvisit_mp(j) == 1
        continue;
    end
    ver_list(end+1) = j;
    isvisit_mp(j) = 1;
    index_vertice = find(labels == labels(j));
    while ( ~isempty(ver_list) ) 
        ver_i = ver_list(1);
        ver_list(1) = [];
        label_mp(ver_i) = label_cnt;
        nei_index = neighbors(:, ver_i);
        nei_index(nei_index < 1) = [];

        index_3 = intersect(nei_index, index_vertice);

        if ~isempty(index_3)
            for k = 1:length(index_3)
                ver_tmp_3 = index_3(k);
                if isvisit_mp(ver_tmp_3) == 0
                    ver_list(end+1) = ver_tmp_3;
                end
                isvisit_mp(ver_tmp_3) = 1;
            end
        end
    end
    label_cnt = label_cnt + 1;
    
end


end