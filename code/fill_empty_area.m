function ind_label = fill_empty_area(ind_label, corr_map, mesh)
% filled the empty area of individual atlas
% Author: Jinlong Li and Guoyuan Yang, BIT.
% 
% Inputs:
%     - ind_label
%       atlas label
% 
%     - corr_map
%       RSFC matrix 
%       
% 
%     - mesh
%       template space. e.g. 'fs_LR_32k'
% 
% Outputs:
%     - ind_label
%       result label

%% prepare
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


index_1 = find(ind_label == 0);
index_1 = intersect(index_1, no_medial_wall_index);

label_mp = zeros(64984,1);
isvisit_mp = ones(64984, 1);
isvisit_mp(index_1) = 0;


label_cnt = 1;
ver_list = [];
    
nei_list_all = {};
nei_list_all{end+1} = [];

for j = 1:length(index_1)

    ver = index_1(j);
    if isvisit_mp(ver) == 1
        continue;
    end
    ver_list(end+1) = ver;
    isvisit_mp(ver) = 1;
    while ( ~isempty(ver_list) ) 
        ver_i = ver_list(1);
        ver_list(1) = [];
        label_mp(ver_i) = label_cnt;
        nei_index = neighbors(:, ver_i);
        nei_index(nei_index < 1) = [];
        index_2 = intersect(nei_index, index_1);
        index_3 = setdiff(nei_index, index_1);
        if ~isempty(index_2)
            for k = 1:length(index_2)
                ver_tmp_2 = index_2(k);
                if isvisit_mp(ver_tmp_2) == 1
                    continue;
                end
                ver_list(end+1) = ver_tmp_2;
                isvisit_mp(ver_tmp_2) = 1;
            end
        end

        if ~isempty(index_3)
            for k = 1:length(index_3)
                ver_tmp_3 = index_3(k);
                label_mp(ver_tmp_3) = 1000+label_cnt;
                nei_list_all{label_cnt} = [nei_list_all{label_cnt}; ind_label(ver_tmp_3)];
            end
        end
    end
    nei_list_all{label_cnt} = unique(nei_list_all{label_cnt});

    label_cnt = label_cnt + 1;
    nei_list_all{end+1} = [];
end



index_2 = find(label_mp > 1000);
a_pa = unique(ind_label(index_2));
a_pa(a_pa == 0) = [];

for j = 1:length(a_pa)
    pa = a_pa(j);
    if j == 1
        index_2 = find(ind_label == pa);
    else
        index_2 = [index_2; find(ind_label == pa)];
    end
end
index_2 = [index_1; index_2 ];



% candidate 100
for j = 1:length(index_1)
    ver_j = index_1(j);
    index_3 = nei_list_all{label_mp(index_1(j))};
    index_3(index_3 < 1) = [];
    index_2 = [];
    for ti = 1:length(index_3)
        index_2 = [index_2; find(ind_label == index_3(ti))];
    end
    
    [~, m_i] = maxk(corr_map(ver_j, index_2), 100);
    ind_label(ver_j) = mode(ind_label(index_2(m_i)));
end

end


