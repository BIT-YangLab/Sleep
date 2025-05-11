function create_MNI_dlabel(label, outfile)
% Create dlabel file in MNI152 with atlas
% Author: Jinlong Li and Guoyuan Yang, BIT.
% 
% Inputs:
%     - label
%       atlas label 
% 
%     - outfile
%       Path to save result

    addpath('/nd_disk3/guoyuan/Xinyu/software/cifti-matlab-master');
    template_dlabel = '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/resource/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_netassignments_LR.dlabel.nii';
    colormap_file = '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/resource/Yeo7_networks_colormap.txt';
    
    template_file_nii = './MNI_immidiate_file.dlabel.nii';
    template_d = cifti_read(template_dlabel);
    load( '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/resource/fs_LR_32k_no_medial_wall_index.mat');  %no_medial_wall_index

    if length(label) == 96854
        label(setdiff(1:64984, no_medial_wall_index)) = [];
    end
    if length(label) ~= 91282
        error('dimension = 91282 or 96854');
    end
    
    template_d.cdata = single(label);
    cifti_write(template_d, template_file_nii);
    system(['wb_command -cifti-label-import ' template_file_nii ' ' colormap_file ' ' outfile]);
    system(['rm ' template_file_nii]);

end