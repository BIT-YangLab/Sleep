function color_label_visual(source, out_dir, file_stem, flg, netw)
% Visualizing cortex
% Author: Jinlong Li and Guoyuan Yang, BIT.
% 
% Inputs:
%     - source
%       Matrix of atlas labels, 59412*1 or 64984*1
% 
%     - out_dir
%       Path to save result
% 
%     - file_stem
%       file stem to save result
% 
%     - flg
%       vertex number
% 
%     - netw
%       number of networks or parcellations
%       
%    
% Outputs:
%     Save to specified file

    load( '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/resource/fs_LR_32k_no_medial_wall_index.mat');  %no_medial_wall_index
    create_label_script='/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/code/Convert_to_cifti_label.py';

    %% color file
%     color_L_file = '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/resource/Lh_Gordon_networks_fs32_new_color.label.gii';
%     color_R_file = '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/resource/Rh_Gordon_networks_fs32_new_color.label.gii';
    if netw == 7
        color_L_file = '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/resource/Lh_7networks_fs32_new.label.gii';
        color_R_file = '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/resource/Rh_7networks_fs32_new.label.gii';
    elseif netw == 14
        color_L_file = '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/resource/Lh_gordon333_13net.label.gii';
        color_R_file = '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/resource/Rh_gordon333_13net.label.gii';
    elseif netw == 17
        color_L_file = '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/resource/Lh_17networks_fs32_new_color.label.gii';
        color_R_file = '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/resource/Rh_17networks_fs32_new_color.label.gii';
    else
        color_L_file = '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/resource/Lh_Gordon_networks_fs32_new_color_fix.label.gii';
        color_R_file = '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/resource/Rh_Gordon_networks_fs32_new_color_fix.label.gii';
    end
 





%     src_label_cifti = ft_read_cifti_mod(source);
%     src_label = src_label_cifti.data;
    if flg == 59412
        new_label = zeros(64984, 1);
        new_label(no_medial_wall_index) = source;
    elseif flg == 64984
        new_label = source;
    end

    lh_labels=new_label(1:32492);
    rh_labels=new_label(32493:64984);
    save([out_dir '/temp_label_1.mat'], 'lh_labels', 'rh_labels');

    L_out_file = [out_dir '/Network_visualization_fill_color_L' file_stem '.label.gii'];
    R_out_file = [out_dir  '/Network_visualization_fill_color_R' file_stem '.label.gii'];

    system(['python  ' create_label_script '   '  out_dir '/temp_label_1.mat  ' color_L_file '  ' color_R_file '  '  L_out_file '  '  R_out_file]);

    
    system(['wb_command -set-structure ' L_out_file ' CORTEX_LEFT']);
    system(['wb_command -set-structure ' R_out_file ' CORTEX_RIGHT']);

end