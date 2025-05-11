% Author: Jinlong Li and Guoyuan Yang, BIT.
% construct individual atlas for subcortex by Lisa2019 in control group(Ji et al., 2015)
% for each subject, construct atlas with single sessions 
clc;clear;

work_dir = '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/';
ret_dir = '/nd_disk3/guoyuan/Jinlong/sleep_ret/HFR_ret/';

% network_label_all sleep_info sublist sub_stage_cell 
load([ret_dir '/IndiPar_net_multi_7/Network_ind_ret5min.mat']);
load([ work_dir '/resource/fs_LR_32k_no_medial_wall_index.mat']);  %no_medial_wall_index
load([ work_dir '/sleep_subinfo.mat']);
% line 6: eo1
sleep_sub_stage_cnt(:, ~filter_sub_idx) = [];
sleep_sub_file_info(:, ~filter_sub_idx) = [];
sleep_sub_stage_info(:, ~filter_sub_idx) = [];
sublist(~filter_sub_idx) = [];
network_label_all(:, ~filter_sub_idx) = [];




load([ ret_dir '/subject_info.mat']);

down_900mesh = ft_read_cifti([ work_dir '/resource/fslr_downsample_900mesh_parcellation.dlabel.nii']);
down_900mesh_label = down_900mesh.dlabel(no_medial_wall_index);

% line 6: eo1

network_n = 7;
vertex_n = 64984;
stageNames = {'Awake', 'N1', 'N2', 'N3', 'REM'};

% Awake N1 N2 N3
network_label_all = network_label_all([6 2 3 4 5], :);


network_label_subcortex_lisa = cell(size(network_label_all));


load([ work_dir '/resource/Template/Yeo' num2str(network_n) 'net_grp_template.mat']);
[mv, mi] = max(group_network_match, [], 1);
group_template_cortex_label = mi - 1;

lh = gifti([ work_dir '/resource/Glasser_2016.32k.L.label.gii']);
rh = gifti([ work_dir '/resource/Glasser_2016.32k.R.label.gii']);
lh = lh.cdata; rh = rh.cdata;
rh = rh + max(lh);
glasser_label = [lh; rh];
glasser_label = glasser_label(no_medial_wall_index);

label_cortex = group_template_cortex_label';
label_cortex = label_cortex(no_medial_wall_index);

glasser_label_list = unique(glasser_label); glasser_label_list(glasser_label_list < 1) = [];
glasser_label_map = zeros(size(glasser_label_list));

for glasser_label_idx = 1:length(glasser_label_list)
    idx = find(glasser_label == glasser_label_list(glasser_label_idx));
    glasser_label_map(glasser_label_idx) = mode(label_cortex(idx));
    
end

load([ work_dir '/resource/sleep_eo_sfnr.mat']);
sfnr_map = sfnr_map(no_medial_wall_index);
sfnr_map = (sfnr_map - min(sfnr_map)) / (max(sfnr_map) - min(sfnr_map));

sfnr_map_glasser = [];
sfnr_map_label = zeros(size(sfnr_map));
for roi_idx = 1:length(glasser_label_list)
    idx = find(glasser_label == glasser_label_list(roi_idx));
    sfnr_map_glasser(roi_idx) = nanmean(sfnr_map(idx));
    sfnr_map_label(idx) = sfnr_map_glasser(roi_idx);
end


Xue_top_x = 400;
Wu_top_x = 20;


for stagei = 1:5
    
    for subi = 1:size(network_label_all, 2)
        disp(['process individualized subcortex maps:  ' num2str(subi) '/' num2str(size(network_label_all, 2))]); tic
        temp_sleep_info = sleep_sub_file_info{stagei, subi};

        if max(size(temp_sleep_info)) == 0

            network_label_subcortex_lisa{stagei, subi} = [];

            continue;
        end

        tmp_label_buckner = {};
        tmp_label_lisa = {};
        tmp_label_wu = {};

        
        for sesi = 1:length(temp_sleep_info)
            cifti_f = temp_sleep_info{sesi};
            raw_BOLD_set = [];
            if ~exist(cifti_f)
                continue;
            end
            cifti_data = ft_read_cifti(cifti_f);
            if stagei == 1
                % For data with eyes open, T1 relaxation gradually reaches a steady state
                % at the beginning of FMRI scanning, and the first few frames are generally
                % not in a stable state. Remove frames around 10 seconds and use the data
                % with steady state afterwards.
                raw_BOLD = cifti_data.dtseries(:, 11:160);
            else
                raw_BOLD = cifti_data.dtseries;
            end
            vol=bsxfun(@minus,raw_BOLD,mean(raw_BOLD,2)); % row wise demeaning
            raw_BOLD=bsxfun(@rdivide,vol,std(vol',1)'); % row wise standardization 
            raw_BOLD_set = [raw_BOLD_set raw_BOLD];
        
    
            if isempty(raw_BOLD_set)
                continue; 
            end

            subcortex_BOLD = raw_BOLD_set(64985:end, :);
            raw_BOLD_set = raw_BOLD_set(no_medial_wall_index, :);
            cortex_BOLD = [];
            for roi_idx = 1:length(glasser_label_list)
                idx = find(glasser_label == glasser_label_list(roi_idx));
                cortex_BOLD(roi_idx, :) = nanmean(raw_BOLD_set(idx, :), 1);
            end

            

            label_cortex = group_template_cortex_label';
            label_cortex = label_cortex(no_medial_wall_index);


            % pearson corr
            corr_mat = paircorr_mod(subcortex_BOLD', cortex_BOLD');
            corr_mat = FisherTransform(corr_mat);
            corr_mat = corr_mat .* (sfnr_map_glasser);
            
            toc;
            
            
            % lisa 2019
            new_corr = [];
            for network_i = 1:network_n
                idx = find(glasser_label_map == network_i);
                new_corr(:, network_i) = mean(corr_mat(:, idx), 2);
            end
            [mv, mi] = max(new_corr, [], 2);
            label_subcortex = (mi);
            label_whole = [network_label_all{stagei, subi}{sesi}; label_subcortex];
            tmp_label_lisa{end+1} = label_whole;

            
        end

        network_label_subcortex_lisa{stagei, subi} = tmp_label_lisa;

        toc;
    end
end