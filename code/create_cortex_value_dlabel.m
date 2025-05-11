function create_cortex_value_dlabel(label, outfile, cmp, min_scale, max_scale)
% create dlabel file for map in cortex
% Author: Jinlong Li and Guoyuan Yang, BIT.
% 
% Inputs:
%     - label
%       map with specific values
% 
%     - outfile
%       Path to save result
%
%     - cmp
%       A n*3 matrix contains RGB color.
%
%     - min_scale/max_scale
%       Range of visualized numerical values
    addpath('/nd_disk3/guoyuan/Xinyu/software/cifti-matlab-master');
    template_dlabel = '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/resource/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_netassignments_LR.dlabel.nii';
    colormap_file_tmp = '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/resource/Temp_networks_colormap.txt';
    
    fid = fopen(colormap_file_tmp, 'w');
    for row_i = 1:size(cmp, 1)
        fwrite(fid, sprintf('%d\n', row_i));
        fwrite(fid, sprintf('%d %d %d %d %d\n', row_i, ceil(cmp(row_i, 1) * 255), ceil(cmp(row_i, 2)*255), ceil(cmp(row_i, 3)*255), 255));
    end
    fclose(fid);
    
    
    template_file_nii = './MNI_immidiate_file.dlabel.nii';
    template_d = cifti_read(template_dlabel);
    load( '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/resource/fs_LR_32k_no_medial_wall_index.mat');  %no_medial_wall_index

    if length(label) == 96854
        label(setdiff(1:64984, no_medial_wall_index)) = [];
    end
    if length(label) ~= 91282
        error('dimension = 91282 or 96854');
    end
    
    new_label = label;
    for row_i = 1:length(label)
        new_label(row_i, :) = ceil((label(row_i)-min_scale)/(max_scale-min_scale) * size(cmp, 1)) + 1;
    end

    template_d.cdata = single(new_label);
    cifti_write(template_d, template_file_nii);
    system(['wb_command -cifti-label-import ' template_file_nii ' ' colormap_file_tmp ' ' outfile]);
    system(['rm ' template_file_nii]);
    system(['rm ' colormap_file_tmp]);
end