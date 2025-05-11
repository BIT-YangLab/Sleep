% Author: Jinlong Li and Guoyuan Yang, BIT.
% compute rsfc between encroaching region for AM network
% 
clc;clear;

work_dir = '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/';
ret_dir = '/nd_disk3/guoyuan/Jinlong/sleep_ret/HFR_ret/';

addpath '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/code'

% network_label_all sleep_info sublist sub_stage_cell 

load([ work_dir '/resource/fs_LR_32k_no_medial_wall_index.mat']);  %no_medial_wall_index
load([ work_dir '/subject_info.mat']);

template_dlabel = '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/resource/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_netassignments_LR.dlabel.nii';

dlabel_struct = ft_read_cifti(template_dlabel);

cortex_mask_lh = dlabel_struct.brainstructure == 1;
cortex_mask_rh = dlabel_struct.brainstructure == 2;
cerebellum_mask_lh = dlabel_struct.brainstructure == 10;
cerebellum_mask_rh = dlabel_struct.brainstructure == 11;
thalamus_mask_lh = dlabel_struct.brainstructure == 20;
thalamus_mask_rh = dlabel_struct.brainstructure == 21;
striatum_mask_lh = dlabel_struct.brainstructure == 3 | dlabel_struct.brainstructure == 8 | dlabel_struct.brainstructure == 18;
striatum_mask_rh = dlabel_struct.brainstructure == 4 | dlabel_struct.brainstructure == 9 | dlabel_struct.brainstructure == 19;


mask_list = {cortex_mask_lh, cortex_mask_rh, cerebellum_mask_lh, cerebellum_mask_rh, thalamus_mask_lh, thalamus_mask_rh, striatum_mask_lh, striatum_mask_rh};
mask_name_list = {'cortex_lh', 'cortex_rh', 'cerebellum_lh', 'cerebellum_rh', 'thalamus_lh', 'thalamus_rh', 'striatum_lh', 'striatum_rh'};
mask_list_w = {cortex_mask_lh | cortex_mask_rh, cerebellum_mask_lh | cerebellum_mask_rh, thalamus_mask_lh | thalamus_mask_rh, striatum_mask_lh | striatum_mask_rh, ...
    cortex_mask_lh | cortex_mask_rh | cerebellum_mask_lh | cerebellum_mask_rh | thalamus_mask_lh | thalamus_mask_rh | striatum_mask_lh | striatum_mask_rh};
mask_name_list_w = {'cortex', 'cerebellum', 'thalamus', 'striatum', 'wholebrain'};

%% select different brain region
mask_idx = 1;

stageNames = {'Awake', 'N1', 'N2', 'N3', 'REM'};

load([ret_dir '/IndiPar_net7/Network_ind_subcortex_cleanup8min.mat']);
network_label_all = network_label_subcortex_lisa;

label_reduction_list = cell(4, 1);

for stagei = 2:4
    label_reduction_stage = zeros(96854, size(network_label_all, 2));
    nc_time = 0;
    for subi = 1:size(network_label_all, 2)
        if isempty(network_label_all{stagei, subi}) || isempty(network_label_all{1, subi}) 
            continue;
        end
        label_1 = network_label_all{stagei, subi}; label_2 = network_label_all{1, subi};
        label_1_idx = find(label_1 == 4); label_2_idx = find(label_2 == 4);
        label_reduction_idx = setdiff(label_2_idx, label_1_idx);
        label_hold_idx = intersect(label_2_idx, label_1_idx);
        label_reduction_stage(label_reduction_idx, subi) = 1;
        label_reduction_stage(label_hold_idx, subi) = 2;
        nc_time = nc_time + 1;
    end
    label_reduction_list{stagei} = label_reduction_stage;
end 

rsfc_stage_list = cell(4, 1);

for stagei = 2:4
    rsfc_tmp_list = zeros(10, 6);
    nc_time = 1;
    for subi = 1:size(sleep_sub_file_info, 2)

        disp(['process individualized subcortex maps:  ' num2str(subi) '/' num2str(size(sleep_sub_file_info, 2))]); tic
        tic;
        if isempty(network_label_all{stagei, subi}) || isempty(network_label_all{1, subi}) 
            continue;
        end

        label_reduction_tmp = label_reduction_list{stagei}(:, subi);
        label_reduction_tmp_list = unique(label_reduction_tmp);
        label_reduction_tmp_list(label_reduction_tmp_list < 1) = [];

        label_reduction_idx = find(label_reduction_tmp ~= 0);
        label_reduction_tmp_new = label_reduction_tmp(label_reduction_idx);

        temp_sleep_info = sleep_sub_file_info{1, subi};
    
        if max(size(temp_sleep_info)) == 0
            continue;
        end
    
        raw_BOLD_set = [];
        for sesi = 1:length(temp_sleep_info)
            cifti_f = temp_sleep_info{sesi};
            
            if ~exist(cifti_f)
                continue;
            end
            cifti_data = ft_read_cifti(cifti_f);
            raw_BOLD = cifti_data.dtseries;
            vol=bsxfun(@minus,raw_BOLD,mean(raw_BOLD,2)); % row wise demeaning
            raw_BOLD=bsxfun(@rdivide,vol,std(vol',1)'); % row wise standardization 
            raw_BOLD_set = [raw_BOLD_set raw_BOLD];
        end
    
        if isempty(raw_BOLD_set)
            continue; 
        end
    
        new_BOLD_set = raw_BOLD_set(label_reduction_idx, :);
        
        corr_mat = paircorr_mod(new_BOLD_set');
        new_corr_mat = zeros(length(label_reduction_tmp_list));
        for ki = 1:length(label_reduction_tmp_list)
            idx_i = find(label_reduction_tmp_new == label_reduction_tmp_list(ki));
            for kj = 1:length(label_reduction_tmp_list)
                idx_j = find(label_reduction_tmp_new == label_reduction_tmp_list(kj));
                new_corr_mat(ki, kj) = mean(mean(corr_mat(idx_i, idx_j)));
            end
        end

        new_corr_mat_awake = new_corr_mat(triu(ones(size(new_corr_mat))) == 1);
    
        
        
        temp_sleep_info = sleep_sub_file_info{stagei, subi};
    
        if max(size(temp_sleep_info)) == 0
            continue;
        end
    
        raw_BOLD_set = [];
        for sesi = 1:length(temp_sleep_info)
            cifti_f = temp_sleep_info{sesi};
            
            if ~exist(cifti_f)
                continue;
            end
            cifti_data = ft_read_cifti(cifti_f);
            raw_BOLD = cifti_data.dtseries;
            vol=bsxfun(@minus,raw_BOLD,mean(raw_BOLD,2)); % row wise demeaning
            raw_BOLD=bsxfun(@rdivide,vol,std(vol',1)'); % row wise standardization 
            raw_BOLD_set = [raw_BOLD_set raw_BOLD];
        end
    
        if isempty(raw_BOLD_set)
            continue; 
        end
    
        new_BOLD_set = raw_BOLD_set(label_reduction_idx, :);
        
        corr_mat = paircorr_mod(new_BOLD_set');
        new_corr_mat = zeros(length(label_reduction_tmp_list));
        for ki = 1:length(label_reduction_tmp_list)
            idx_i = find(label_reduction_tmp_new == label_reduction_tmp_list(ki));
            for kj = 1:length(label_reduction_tmp_list)
                idx_j = find(label_reduction_tmp_new == label_reduction_tmp_list(kj));
                new_corr_mat(ki, kj) = mean(mean(corr_mat(idx_i, idx_j)));
            end
        end
        new_corr_mat = new_corr_mat(triu(ones(size(new_corr_mat))) == 1);
        
        rsfc_tmp_list(nc_time, :) = [new_corr_mat_awake; new_corr_mat];
        nc_time = nc_time + 1;
        toc;
    end
    rsfc_stage_list{stagei} = rsfc_tmp_list;
end


