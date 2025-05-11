% Author: Jinlong Li and Guoyuan Yang, BIT.
% visualization for group average atlas 
% cerebral cortex: dlabel.nii, cerebellum: suit, striatum&thalamus: scatter point
clc;clear;

work_dir = '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/';
ret_dir = '/nd_disk3/guoyuan/Jinlong/sleep_ret/HFR_ret/';

addpath '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/code'

% network_label_all sleep_info sublist sub_stage_cell 

load([ work_dir '/resource/fs_LR_32k_no_medial_wall_index.mat']);  %no_medial_wall_index
load([ work_dir '/sleep_subinfo.mat']);



stageNames = {'Awake', 'N1', 'N2', 'N3', 'REM'};

load([ret_dir '/IndiPar_net7/Network_ind_subcortex_cleanup8min.mat']);
% Yeo_7Networks_Color
load([work_dir '/resource/Yeo_7Networks_Color.mat']);

%% clean up
template_dlabel = [ work_dir '/resource/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_netassignments_LR.dlabel.nii'];
template_MNI152_cortex = [ work_dir '/resource/MNI152_cortex_mask.nii.gz'];

a = ft_read_cifti(template_dlabel);
cortex_mask_lh = a.brainstructure == 1;
cortex_mask_rh = a.brainstructure == 2;
cerebellum_mask_lh = a.brainstructure == 10;
cerebellum_mask_rh = a.brainstructure == 11;
thalamus_mask_lh = a.brainstructure == 20;
thalamus_mask_rh = a.brainstructure == 21;
striatum_mask_lh = a.brainstructure == 3 | a.brainstructure == 8 | a.brainstructure == 18;
striatum_mask_rh = a.brainstructure == 4 | a.brainstructure == 9 | a.brainstructure == 19;

mask_list = {cortex_mask_lh, cortex_mask_rh, cerebellum_mask_lh, cerebellum_mask_rh, thalamus_mask_lh, thalamus_mask_rh, striatum_mask_lh, striatum_mask_rh};
mask_list_w = {cortex_mask_lh | cortex_mask_rh, cerebellum_mask_lh | cerebellum_mask_rh, thalamus_mask_lh | thalamus_mask_rh, striatum_mask_lh | striatum_mask_rh, ...
    cortex_mask_lh | cortex_mask_rh | cerebellum_mask_lh | cerebellum_mask_rh | thalamus_mask_lh | thalamus_mask_rh | striatum_mask_lh | striatum_mask_rh};
mask_name_list_w = {'cortex', 'cerebellum', 'thalamus', 'striatum', 'wholebrain'};


group_label_list = cell(5, 1);



subcortex_label = network_label_subcortex_lisa;
    

group_network_size_list = zeros(length(mask_list), 4);
for stagei = 1:5
    label1 = [];
    for si = 1:size(subcortex_label, 2)
        if isempty(subcortex_label{stagei, si})
            continue;
        end
    
        label1 = [label1 subcortex_label{stagei, si}];
    end
    label_average = mode(label1, 2);
    group_label_list{stagei} = label_average;

    
    plot_subcortex(label_average, [ret_dir '/buckner_test.nii.gz'], [ret_dir '/pic/Subcortex_' subcortex_name '_stage_' stageNames{stagei}]);
    create_MNI_dlabel(label_average, [ret_dir '/MNI152_subcortex_' subcortex_name '_visualization_stage_' stageNames{stagei} '.dlabel.nii']);
    create_cortex_label_dlabel(label_average, [ret_dir '/sleep_network_' stageNames{stagei}  '.dlabel.nii'], Yeo_7Networks_Color);
end



%% hemisphere group network size
group_size_list = zeros(4, 10);
group_size_list_w = zeros(4, 5);
for mask_idx = 1:4
    
    for stagei = 1:5
        label_tmp_lh = group_label_list{stagei}(mask_list{1+(mask_idx-1)*2});
        label_tmp_rh = group_label_list{stagei}(mask_list{mask_idx*2});
        label_tmp_all = group_label_list{stagei}(mask_list_w{mask_idx});
        group_size_list(mask_idx, 1+(stagei-1) * 2) = nnz(label_tmp_lh==4) / length(label_tmp_lh);
        group_size_list(mask_idx, stagei * 2) = nnz(label_tmp_rh==4) / length(label_tmp_rh);
        group_size_list_w(mask_idx, stagei) = nnz(label_tmp_all == 4) / length(label_tmp_all);
    end
    
end

% for better visualization, we remove the vertexes of medial wall when
% computing group network size in cortex.
network_all_vertex = [59412, nnz(mask_list_w{2}), nnz(mask_list_w{3}), nnz(mask_list_w{4}), 59412+nnz(mask_list_w{2})+nnz(mask_list_w{3})+nnz(mask_list_w{4})];

group_size_all_network_w = cell(5, 1);
for mask_idx = 1:5
    group_size_all_network_w_tmp = zeros(4, 7);
    for networki = 1:7
        for stagei = 1:4
            
            label_tmp_all = group_label_list{stagei}(mask_list_w{mask_idx});
            group_size_all_network_w_tmp(stagei, networki) = nnz(label_tmp_all == networki) / network_all_vertex(mask_idx);
        end
    end
    group_size_all_network_w{mask_idx} = group_size_all_network_w_tmp;
end

group_size_list_v = zeros(4, 4);
for i = 1:4
    group_size_list_v(:, i) = (group_size_list_w(:, i+1) - group_size_list_w(:, 1)) ./  group_size_list_w(:, 1);
end


