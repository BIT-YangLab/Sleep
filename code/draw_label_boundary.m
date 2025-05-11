function labels2 = draw_label_boundary(new_label, mesh)
% find boundary between networks or parcellations in atlas
% Author: Jinlong Li and Guoyuan Yang, BIT.
% 
% Inputs:
%     - new_label
%       atlas lable
% 
%     - mesh
%       template space. e.g. 'fs_LR_32k'
% Outputs:
%     - labels2
%       result label


%% labels2 = draw_label_boundary(new_label)
% new_label -> vertices * 1 : label matrix
% mesh: example: fs_LR_32k
CBIG_PATH = '/nd_disk3/guoyuan/Jinlong/CBIG/';
addpath([CBIG_PATH '/utilities/matlab/fslr_matlab/']);
addpath([CBIG_PATH '/external_packages/SD/SDv1.5.1-svn593/BasicTools/']);


if contains(mesh, 'fs_LR')
    lh_targ_mesh = CBIG_read_fslr_surface('lh', 'fs_LR_32k','inflated','medialwall.annot');
    rh_targ_mesh = CBIG_read_fslr_surface('rh', 'fs_LR_32k','inflated','medialwall.annot');
    no_medial_wall_index = [find(lh_targ_mesh.MARS_label == 2); find(rh_targ_mesh.MARS_label == 2) + length(lh_targ_mesh.MARS_label)];
    lh_end_idx = 32492;
elseif contains(mesh, 'fsaverage') 
    lh_targ_mesh = CBIG_ReadNCAvgMesh('lh', mesh, 'inflated', 'cortex');
    rh_targ_mesh = CBIG_ReadNCAvgMesh('rh', mesh, 'inflated', 'cortex');
    
    no_medial_wall_index = [find(lh_targ_mesh.MARS_label == 2) find(rh_targ_mesh.MARS_label == 2) + length(lh_targ_mesh.MARS_label)]';
    lh_end_idx = 40962;
end

new_label1 = new_label;

lh_label = new_label1(1:lh_end_idx);
rh_label = new_label1(lh_end_idx+1:end);

lh_list = unique(lh_label);
lh_list(lh_list < 1) = [];
rh_list = unique(rh_label);
rh_list(rh_list < 1) = [];

for pa = 1:length(lh_list)
    index = find(lh_label == lh_list(pa));
    for i = 1:length(index)
        neighbor_index = lh_targ_mesh.vertexNbors(:, index(i));
        neighbor_index(neighbor_index < 1) = [];
        if isempty(neighbor_index)
            continue
        end
        el_index = intersect(setdiff(neighbor_index, index), no_medial_wall_index);
        tmp_label = lh_label(el_index);
        if ~isempty(el_index) && nnz(tmp_label) ~= 0
            lh_label(index(i)) = 0;
        end 
    end
end

for pa = 1:length(rh_list)
    index = find(rh_label == rh_list(pa));
    for i = 1:length(index)
        neighbor_index = rh_targ_mesh.vertexNbors(:, index(i));
        neighbor_index(neighbor_index < 1) = [];
        if isempty(neighbor_index)
            continue
        end
        el_index = intersect(setdiff(neighbor_index, index), no_medial_wall_index-lh_end_idx);
        tmp_label = rh_label(el_index);
        if ~isempty(el_index) && nnz(tmp_label) ~= 0
            rh_label(index(i)) = 0;
        end 
    end
end

labels2 = [lh_label; rh_label];

end