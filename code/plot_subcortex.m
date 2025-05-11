function plot_subcortex(label, nii_file, out_fig_base_dir)
% plot subcortex
% Author: Xinyu Wu and Guoyuan Yang, BIT.
% 
% 
    addpath(genpath('/nd_disk3/guoyuan/Xinyu/plot_fig_subcortex'));
    addpath('/nd_disk3/guoyuan/Xinyu/software/cifti-matlab-master');
    addpath('/nd_disk3/guoyuan/Xinyu/software/spm12');

    %% Convert results into nifti files.
    cifti_base_info = ft_read_cifti('/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/resource/CortexSubcortex_ColeAnticevic_NetPartition_woSubcorGSR_netassignments_LR.dscalar.nii');
    nifti_base_info = MRIread('/nd_disk3/guoyuan/Xinyu/plot_fig_subcortex/rois/Atlas_ROIs.2.nii.gz');
    x = 63 + cifti_base_info.pos(64985:96854, 2)/2;
    y = 45 - cifti_base_info.pos(64985:96854, 1)/2;
    z = 36 + cifti_base_info.pos(64985:96854, 3)/2;
    pos = sub2ind([109, 91, 91], x, y, z);

    subcortex_label = label(64985:end);
    data = zeros(109, 91, 91);
    data(pos) = subcortex_label;
    nifti_base_info.vol = data;
    MRIwrite(nifti_base_info, nii_file);

    min_scale = 0;
    max_scale = 7;
    load('/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/resource/Yeo_7Networks_Color.mat');
    cmap = Yeo_7Networks_Color;%flip(slanCM(100, 256));

     age_grad_hiptha_png_path = [out_fig_base_dir '_hiptha.png'];
     plot_fig_subcortex_mod(nii_file, 'export_fig_path', age_grad_hiptha_png_path, 'cmap', cmap, 'cscale', [min_scale, max_scale], 'nifti_roi_label', [10, 49]);
     camorbit(191, 4, 'data', [0, 0, 1]);
%     camorbit(95, -8, 'data', [0, 0, 1]);
     export_fig(age_grad_hiptha_png_path, '-m6', '-q100');
     close;
            
    age_grad_str_png_path = [out_fig_base_dir '_str.png'];
    plot_fig_subcortex_mod(nii_file, 'export_fig_path', age_grad_str_png_path, 'cmap', cmap, 'cscale', [min_scale, max_scale], 'nifti_roi_label', [11, 12, 13, 26, 50, 51, 52, 58]);
    camorbit(39, 0, 'data', [0, 0, 1]);
    export_fig(age_grad_str_png_path, '-m6', '-q100');
    close;

    age_grad_cbm_png_path = [out_fig_base_dir '_cbm.png'];
    map = suit_map2surf(nii_file, 'space', 'SPM', 'stats',@mode);
    set(gcf, 'Position', [0, 0, 1000, 1000]);
    set(gca, 'xticklabel', []); set(gca, 'yticklabel', []); set(gca, 'zticklabel', []);
    set(gca, 'color', 'none'); set(gcf, 'color', 'none');
    axis off;
    suit_plotflatmap(map, 'cmap', cmap, 'cscale', [min_scale, max_scale], 'type','label');
    export_fig(age_grad_cbm_png_path, '-m6', '-q100');
    close;

end