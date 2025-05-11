% Author: Jinlong Li and Guoyuan Yang, BIT.
% compute sfnr map for subcortex mapping
% 
clc;clear;

work_dir = '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/';
ret_dir = '/nd_disk3/guoyuan/Jinlong/sleep_ret/HFR_ret/';

% network_label_all sleep_info sublist sub_stage_cell 
load([ret_dir '/IndiPar_net7/Network_ind_ret8min.mat']);
load([ work_dir '/resource/fs_LR_32k_no_medial_wall_index.mat']);  %no_medial_wall_index
load([ work_dir '/sleep_subinfo.mat']);
% line 6: eo1
sleep_sub_stage_cnt(:, 34:66) = [];
sleep_sub_file_info(:, 34:66) = [];
sleep_sub_stage_info(:, 34:66) = [];
sublist(34:66) = [];
network_label_all(:, 34:66) = [];



cifti_format = '/nd_disk3/guoyuan/sleep/sleep_eo_ciftify_%s/out/ciftify/sub-%d/MNINonLinear/Results/ses-01_task-sleep_desc-preproc/ses-01_task-sleep_desc-preproc_Atlas_s0.dtseries.nii';
% 
for stagei = 1:size(sleep_sub_file_info, 1)
    for subi = 1:size(sleep_sub_file_info, 2)
        if isempty(sleep_sub_file_info{stagei, subi})
            continue;
        end
        tmp_file_list = {};
        for sesi = 1:length(sleep_sub_file_info{stagei, subi})
            file1 = sleep_sub_file_info{stagei, subi}{sesi};
            cifti_f1 = strsplit(file1, '2.55');
            cifti_f = [cifti_f1{1} '8mm_subcortex_regress.dtseries.nii'];
            if ~exist(cifti_f)
                disp(cifti_f);
            end
            tmp_file_list{end+1} = cifti_f;
        end
        %sleep_sub_file_info{stagei, subi} = tmp_file_list;
    end
end
sfnr_map = [];
nc = 0;
for stagei = 6:6
    for subi = 1:size(sleep_sub_file_info, 2)
        disp(['process ' num2str(subi) '/' num2str(size(sleep_sub_file_info, 2))]);tic;
        if isempty(sleep_sub_file_info{stagei, subi})
            continue;
        end
        for sesi = 1:length(sleep_sub_file_info{stagei, subi})
            tmp_sess = sleep_sub_stage_info{stagei, subi}{sesi};
            tmp_file = sprintf(cifti_format, tmp_sess, sublist(subi));
            if exist(tmp_file)
               
                cifti_s = ft_read_cifti(tmp_file);
                BOLD_data = cifti_s.dtseries;
                sfnr1 = compute_sfnr(BOLD_data);
                if isempty(sfnr_map)
                    sfnr_map = sfnr1;
                else
                    sfnr_map = sfnr_map + sfnr1;
                end
                nc = nc + 1;
            end
        end
        toc;
    end
end
sfnr_map = sfnr_map ./ nc;
save('/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/resource/sleep_eo_sfnr.mat', 'sfnr_map');