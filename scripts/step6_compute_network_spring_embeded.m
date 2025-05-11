clc;clear;

work_dir = '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/';
ret_dir = '/nd_disk3/guoyuan/Jinlong/sleep_ret/HFR_ret';

load('/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/subject_info.mat');
load([ work_dir '/resource/fs_LR_32k_no_medial_wall_index.mat']);  %no_medial_wall_index

stage_names = {'Awake', 'N1', 'N2', 'N3', 'REM'};

%% glasser 360
lh = gifti([ work_dir '/resource/Glasser_2016.32k.L.label.gii']);
rh = gifti([ work_dir '/resource/Glasser_2016.32k.R.label.gii']);
lh = lh.cdata; rh = rh.cdata;
rh = rh + max(lh);
glasser_label = [lh; rh];
glasser_label = glasser_label(no_medial_wall_index);
glasser_label_list = unique(glasser_label); glasser_label_list(glasser_label_list < 1) = [];


%% cluster subcortex
template_dlabel = '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/resource/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_netassignments_LR.dlabel.nii';
template_MNI152_cortex = '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/resource/MNI152_cortex_mask.nii.gz';

% dist_subcortex_subcortex
dlabel_struct = ft_read_cifti(template_dlabel);
% dist_subcortex_subcortex = compute_vertices_dist(dlabel_struct.pos(64985:end, :), dlabel_struct.pos(64985:end, :));
% dist_matrix = zeros(96854, 96854);
% dist_matrix(64985:end, 64985:end) = dist_subcortex_subcortex;

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
mask_list_w = {cortex_mask_lh | cortex_mask_rh,  ...
     cerebellum_mask_lh | cerebellum_mask_rh | thalamus_mask_lh | thalamus_mask_rh | striatum_mask_lh | striatum_mask_rh};
mask_name_list_w = {'cortex', 'wholesubcortex'};





load('/nd_disk3/guoyuan/Jinlong/sleep_ret/HFR_ret/IndiPar_net7/group_atlas.mat');
mid_L = '/nd_disk3/guoyuan/Jinlong/Gordon2016_2/Parcels/Conte69.L.midthickness.32k_fs_LR.surf.gii';
mid_R = '/nd_disk3/guoyuan/Jinlong/Gordon2016_2/Parcels/Conte69.R.midthickness.32k_fs_LR.surf.gii';
dt = ft_read_cifti('/nd_disk3/guoyuan/Jinlong/test.dtseries.nii');
p = zeros(7, 1);



flg_exist_subcortex_roi = 0;

if flg_exist_subcortex_roi == 0
    cluster_label = cell(2, 1);
    match_label_list = cell(2, 1);
    for mask_idx = 1:2
        new_label = zeros(96854, 7);
        for ki = 1:7
            ak = group_label_list{1};
            ak(ak ~= ki) = 0;
            ak(~mask_list_w{mask_idx}) = 0;
            dt.dtseries(:, 1) = ak;
            ft_write_cifti('/nd_disk3/guoyuan/Jinlong/test', dt, 'parameter', 'dtseries');
            system(['wb_command -cifti-find-clusters '  '/nd_disk3/guoyuan/Jinlong/test.dtseries.nii  0 0 0 0 COLUMN '  '/nd_disk3/guoyuan/Jinlong/test_ret.dtseries.nii -left-surface ' ...
                mid_L ' -right-surface ' mid_R ]);
            cluster_dt = ft_read_cifti(['/nd_disk3/guoyuan/Jinlong/test_ret.dtseries.nii']);
            c1 = cluster_dt.dtseries(mask_list_w{mask_idx});
            c2 = unique(c1);
            for bi = 1:length(c2)
                idx = find(c1 == c2(bi));
                if length(idx) < 4
                    c1(idx) = 0;
                end
            end
            c2 = unique(c1);
            c2(c2 < 1) = [];
            zp(ki) = length(c2);
            new_label(mask_list_w{mask_idx}, ki) = c1;
        end
        n1 = zeros(96854, 1);
        for ki = 1:7
            idx = find(new_label(:, ki));
            start_idx = max(n1);
            n1(idx) = start_idx + new_label(idx, ki);
        end

        label_l = unique(n1); label_l(label_l < 1 ) = [];
        group_cluster_1 = zeros(96854, 1);
        roi_network_match1 = zeros(size(label_l));
        for ki = 1:length(label_l)
            idx = find(n1 == label_l(ki));
            group_cluster_1(idx) = ki;
            roi_network_match1(ki) = mode(group_label_list{1}(idx));
        end

        cluster_label{mask_idx} = group_cluster_1;
        match_label_list{mask_idx} = roi_network_match1;
    end

    cluster_label_all = cluster_label{1};
    cluster_label_subcortex = max(cluster_label_all) + cluster_label{2};
    cluster_label_subcortex(cluster_label{2} == 0) = 0;
    cluster_label_all(64985:end) = cluster_label_subcortex(64985:end);

    match_label_list_all = [match_label_list{1}; match_label_list{2}];

    cluster_label_all_list = unique(cluster_label_all); cluster_label_all_list(cluster_label_all_list < 1) = [];
    save([ret_dir '/RSFC_cluster_subcortex_grp_label.mat'], 'cluster_label_subcortex');
else
    load([ret_dir '/RSFC_cluster_subcortex_grp_label.mat']);
end


subcortex_cluster_label = cluster_label_subcortex;
subcortex_cluster_label = subcortex_cluster_label(64985:end);
subcortex_cluster_label_list = unique(subcortex_cluster_label); subcortex_cluster_label_list(subcortex_cluster_label_list < 1) = [];

% group average atlas for Awake
load([ ret_dir '/IndiPar_net7/group_atlas.mat']);

group_average_atlas = group_label_list{1};

group_average_cortex = group_average_atlas(no_medial_wall_index);
group_average_subcortex = group_average_atlas(64985:96854);

% AMN
network_specific = 4;
network_n = 7;


cortex_RSFC_roi_cell = cell(5, 1);
cortex_RSFC_network_cell = cell(5, 1);
subcortex_RSFC_roi_cell = cell(5, 1);
subcortex_RSFC_network_cell = cell(5, 1);

inter_RSFC_roi_cell = cell(5, 1);

match_list_cortex = zeros(length(glasser_label_list), 1);
for ki = 1:length(glasser_label_list)
    idx = find(glasser_label == glasser_label_list(ki));
    match_list_cortex(ki) = mode(group_average_cortex(idx));
end

match_list_subcortex = zeros(length(subcortex_cluster_label_list), 1);
for ki = 1:length(subcortex_cluster_label_list)
    idx = find(subcortex_cluster_label == subcortex_cluster_label_list(ki));
    match_list_subcortex(ki) = mode(group_average_subcortex(idx));
end

flg_exist_rsfc = 0;

if flg_exist_rsfc == 0

for stagei = 2:2

    cortex_RSFC_roi_roi = cell(size(sleep_sub_file_info, 2), 1);
    cortex_RSFC_roi_network = zeros(length(glasser_label_list), network_n);
    subcortex_RSFC_roi_roi = cell(size(sleep_sub_file_info, 2), 1);
    subcortex_RSFC_roi_network = zeros(length(subcortex_cluster_label_list), network_n);
    inter_RSFC_roi_roi = cell(size(sleep_sub_file_info, 2), 1);
    nc = 0;

    for subi = 1:size(sleep_sub_file_info, 2)
        disp(['process stage:' num2str(stagei) ' sub:' num2str(subi)]); tic;
        if isempty(sleep_sub_file_info{stagei, subi})
            continue;
        end

        temp_sleep_info = sleep_sub_file_info{stagei, subi};

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

        new_BOLD_set_cortex = raw_BOLD_set(no_medial_wall_index, :);
        new_BOLD_set_subcortex = raw_BOLD_set(64985:96854, :);

        % pearson correlation and Z-transform
        rsfc_cortex = paircorr_mod(new_BOLD_set_cortex');
        rsfc_cortex(isnan(rsfc_cortex)) = 0;
        rsfc_cortex = single(FisherTransform(rsfc_cortex));

        rsfc_subcortex = paircorr_mod(new_BOLD_set_subcortex');
        rsfc_subcortex(isnan(rsfc_subcortex)) = 0;
        rsfc_subcortex = single(FisherTransform(rsfc_subcortex));

        rsfc_cortex_subcortex = paircorr_mod(new_BOLD_set_cortex', new_BOLD_set_subcortex');



        new_fc_cortex_roi = zeros(length(glasser_label_list));
        new_fc_cortex_network = zeros(length(glasser_label_list), network_n);

        new_fc_subcortex_roi = zeros(length(subcortex_cluster_label_list));
        new_fc_subcortex_network = zeros(length(subcortex_cluster_label_list), network_n);

        new_fc_cortex_subcortex_roi = zeros(length(glasser_label_list), length(subcortex_cluster_label_list));


        % glasser roi-roi
        for li = 1:length(glasser_label_list)
            idx_i = find(glasser_label == glasser_label_list(li));
            for lj = li+1:length(glasser_label_list)
                idx_j = find(glasser_label == glasser_label_list(lj));
                new_fc_cortex_roi(li, lj) = nanmean(nanmean(rsfc_cortex(idx_i, idx_j)));
            end
            for lj = 1:network_n
                idx_j = find(group_average_cortex == lj);
                new_fc_cortex_network(li, lj) = nanmean(nanmean(rsfc_cortex(idx_i, idx_j)));
            end
        end

        for li = 1:length(subcortex_cluster_label_list)
            idx_i = find(subcortex_cluster_label == subcortex_cluster_label_list(li));
            for lj = li+1:length(subcortex_cluster_label_list)
                idx_j = find(subcortex_cluster_label == subcortex_cluster_label_list(lj));
                new_fc_subcortex_roi(li, lj) = nanmean(nanmean(rsfc_subcortex(idx_i, idx_j)));
            end
            for lj = 1:network_n
                idx_j = find(group_average_subcortex == lj);
                new_fc_subcortex_network(li, lj) = nanmean(nanmean(rsfc_subcortex(idx_i, idx_j)));
            end
        end

        for li = 1:length(glasser_label_list)
            idx_i = find(glasser_label == glasser_label_list(li));
            for lj = 1:length(subcortex_cluster_label_list)
                idx_j = find(subcortex_cluster_label == subcortex_cluster_label_list(lj));
                new_fc_cortex_subcortex_roi(li, lj) = nanmean(nanmean(rsfc_cortex_subcortex(idx_i, idx_j)));
            end
            
        end

        nc = nc + 1;
         cortex_RSFC_roi_roi{subi} = new_fc_cortex_roi;
%         cortex_RSFC_roi_network = cortex_RSFC_roi_network + new_fc_cortex_network;
        subcortex_RSFC_roi_roi{subi} = new_fc_subcortex_roi;
%         subcortex_RSFC_roi_network = subcortex_RSFC_roi_network + new_fc_subcortex_network; 
        inter_RSFC_roi_roi{subi} = new_fc_cortex_subcortex_roi;
        toc;
    end
     cortex_RSFC_roi_cell{stagei} = cortex_RSFC_roi_roi ;
%     cortex_RSFC_network_cell{stagei} = cortex_RSFC_roi_network ;
    subcortex_RSFC_roi_cell{stagei} = subcortex_RSFC_roi_roi ;
%     subcortex_RSFC_network_cell{stagei} = subcortex_RSFC_roi_network ;
    inter_RSFC_roi_cell{stagei} = inter_RSFC_roi_roi;
    save([ret_dir '/RSFC_spring_embeded_cluster_stage' num2str(stagei) '.mat'], 'subcortex_RSFC_roi_cell', 'inter_RSFC_roi_cell');
end
else
    load([ret_dir '/RSFC_spring_embeded.mat']);
end

RSFC_roi_average_cell = cell(5, 1);
for stagei = 1:5
    RSFC_roi_average_tmp = zeros(length(glasser_label_list)+length(subcortex_cluster_label_list));
    nc = 0;
    for subi = 1:size(cortex_RSFC_roi_cell{stagei}, 1)
        if isempty(cortex_RSFC_roi_cell{stagei}{subi})
            continue;
        end
        idx_cortex = 1:length(glasser_label_list);
        idx_subcortex = (length(glasser_label_list)+1) : (length(glasser_label_list) + length(subcortex_cluster_label_list));
        RSFC_roi_average_tmp(idx_cortex, idx_cortex) = RSFC_roi_average_tmp(idx_cortex, idx_cortex) + ...
            cortex_RSFC_roi_cell{stagei}{subi}(idx_cortex, idx_cortex);
        RSFC_roi_average_tmp(idx_cortex, idx_subcortex) = RSFC_roi_average_tmp(idx_cortex, idx_subcortex) + ...
            cortex_RSFC_roi_cell{stagei}{subi}(idx_cortex, idx_subcortex);
        RSFC_roi_average_tmp(idx_subcortex, idx_subcortex) = RSFC_roi_average_tmp(idx_subcortex, idx_subcortex) + ...
            cortex_RSFC_roi_cell{stagei}{subi}(idx_subcortex, idx_subcortex);
        nc = nc + 1;
    end
    RSFC_roi_average_tmp = RSFC_roi_average_tmp ./ nc;
    RSFC_roi_average_cell{stagei} = RSFC_roi_average_tmp;
end

% RSFC_roi_average_cell{stagei}: 5x1 cell, the group average of RSFC_roi_cell 
for stagei = 1:5
    cortex_RSFC_roi_roi = RSFC_roi_average_cell{stagei}(1:360, 1:360);
    cortex_RSFC_roi_roi(eye(size(cortex_RSFC_roi_roi, 1)) == 1) = 0;
    match_list_cortex = match_list_r2n(1:360);

    subcortex_RSFC_roi_roi = RSFC_roi_average_cell{stagei}(361:end, 361:end);
    subcortex_RSFC_roi_roi(eye(size(subcortex_RSFC_roi_roi, 1)) == 1) = 0;
    match_list_subcortex = match_list_r2n(361:end);

    all_RSFC_roi_roi = RSFC_roi_average_cell{stagei};
    all_RSFC_roi_roi(eye(size(all_RSFC_roi_roi, 1) == 1)) = 0;
    match_list_all = match_list_r2n;
    match_list_all(361:end) = match_list_all(361:end) + 7;

    edge_list = zeros(1, 3);
    direc_type = cell(1, 1);
    nc = 0;
    for roi_i = 1:size(all_RSFC_roi_roi, 1)
        
        for roi_j = roi_i+1:size(all_RSFC_roi_roi, 2)
            
                nc = nc + 1;
                edge_list(nc, :) = [roi_i, roi_j, all_RSFC_roi_roi(roi_i, roi_j)]; 
                direc_type{nc} = 'Undirected';
            
        end
    end


    edge_list = table(edge_list(:, 1), edge_list(:, 2), direc_type', edge_list(:, 3), 'VariableNames', ...
    {'Source', 'Target', 'Type', 'Weight'});

    [ia, ib] = sort(edge_list.Weight, 'descend');
    ths = ia(ceil(length(ia)*threshold_conn));
    edge_list(edge_list.Weight < ths, :) = []; 
    writetable(edge_list, [ret_dir '/RSFC_cluster_all_edge_list_ths_' num2str(threshold_conn*100) '_stage_' stage_names{stagei} '.csv']);
    
    node_list = table([1:length(match_list_all)]', match_list_all, 'VariableNames', {'Id', 'Network'});
    writetable(node_list, [ret_dir '/RSFC_cluster_all_node_list_ths_' num2str(threshold_conn*100) 'stage_' stage_names{stagei} '.csv']);



end