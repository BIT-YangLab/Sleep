clc;clear;

work_dir = '/home/jinlong/VDisk1/Jinlong/2025_HFR_ai_mod/';
ret_dir = '/home/jinlong/VDisk1/Jinlong/2025_HFR_ai_mod/result';

load([work_dir '/subject_info.mat']);
load([ work_dir '/resource/fs_LR_32k_no_medial_wall_index.mat']);  %no_medial_wall_index

stage_names = {'Awake', 'N1', 'N2', 'N3', 'REM'};

template_dlabel = [ work_dir '/resource/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_netassignments_LR.dlabel.nii'];
template_MNI152_cortex = [ work_dir '/resource/MNI152_cortex_mask.nii.gz' ];

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
mask_list_w = {cortex_mask_lh | cortex_mask_rh, cerebellum_mask_lh | cerebellum_mask_rh,  ...
     cerebellum_mask_lh | cerebellum_mask_rh | thalamus_mask_lh | thalamus_mask_rh | striatum_mask_lh | striatum_mask_rh};
mask_name_list_w = {'cortex', 'cerebellum', 'whole_subcortex'};

% group average atlas for Awake
load([ ret_dir '/IndiPar_net7/group_atlas.mat']);

group_average_atlas = group_label_list{1};

group_average_cortex = group_average_atlas(no_medial_wall_index);
group_average_subcortex = group_average_atlas(64985:96854);

% glasser 360
lh = gifti([ work_dir '/resource/Glasser_2016.32k.L.label.gii']);
rh = gifti([ work_dir '/resource/Glasser_2016.32k.R.label.gii']);
lh = lh.cdata; rh = rh.cdata;
rh = rh + max(lh);
glasser_label = [lh; rh];
glasser_label = glasser_label(no_medial_wall_index);
glasser_label_list = unique(glasser_label); glasser_label_list(glasser_label_list < 1) = [];

% tian 2022 scale-IV
load([ret_dir '/RSFC_cluster_subcortex_grp_label.mat']);
subcortex_label = cluster_label_subcortex(64985:end);
subcortex_label_list = unique(subcortex_label); subcortex_label_list(subcortex_label_list < 1) = [];

subcortex_label_roi_filter = [];
for ki = 1:length(subcortex_label_list)
    idx = find(subcortex_label == subcortex_label_list(ki));
    if nnz(mask_list_w{3}(idx + 64984))
        subcortex_label_roi_filter(end+1) = subcortex_label_list(ki);
    end
end

% AMN
network_specific = 4;
network_n = 7;




load([ret_dir '/RSFC_spring_embeded_cluster.mat']);
load([ret_dir '/IndiPar_net7/Network_ind_subcortex_cleanup8min.mat']);

load([work_dir '/subject_info.mat']);

AgeInYears = Ancillary_data.AgeInYears;
Gender = Ancillary_data.Gender;
EducationYear = Ancillary_data.EducationInYears;



%% individual atlas all
kidx = [1:360 360+subcortex_label_roi_filter];
M_list_cortex_cell = cell(size(network_label_subcortex_lisa, 2), 1);
M_all_list = zeros(size(network_label_subcortex_lisa, 2), 5);

for subi = 1:length(RSFC_roi_cell{stagei})
    M_list = zeros(length(0.15:0.01:0.15), 5);
    if isempty(network_label_subcortex_lisa{1, subi})
        continue;
    end
    label_ind = network_label_subcortex_lisa{1, subi};
    label_ind_cortex = label_ind(no_medial_wall_index);
    label_ind_subcortex = label_ind(64985:end);
    match_list_cortex = zeros(length(glasser_label_list), 1);
    for ki = 1:length(glasser_label_list)
        idx = find(glasser_label == glasser_label_list(ki));
        match_list_cortex(ki) = mode(label_ind_cortex(idx));
    end
    match_list_subcortex = zeros(length(subcortex_label_roi_filter), 1);
    for ki = 1:length(subcortex_label_roi_filter)
        idx = find(subcortex_label == subcortex_label_roi_filter(ki));
        match_list_subcortex(ki) = mode(label_ind_subcortex(idx)) ;
    end
    new_l = [match_list_cortex; match_list_subcortex]; new_l(new_l ~= 4) = 11;

    for stagei = 1:5
    
        if isempty(RSFC_roi_cell{stagei}{subi})
            continue;
        end

        r1 = RSFC_roi_cell{stagei}{subi}(kidx, kidx);
        nc = 0;
        for c1 = 0.15:0.01:0.15
            nc = nc + 1;
            r1_r = r1;
            [ia, ib] = sort(r1_r(:), 'descend'); ths = ia(ceil(length(ia)*c1));
            r1_r(r1_r < ths) = 0; r1_r(r1_r >= ths) = 1;
            M = compute_modularity(r1_r, new_l); M_list(nc, stagei) = M;
        end
    end
    M_list_cortex_cell{subi} = M_list;
    M_all_list(subi, :) = mean(M_list, 1);
end


%% LME create 

variable_Names_list = {'NS', 'Stage', 'Age', 'Gender', 'Education', 'participant_id', 'epoch', 'fd'};



Modularity_all_cell = [];
nc_time = 0;
for stagei = 1:5
    for subi = 1:size(sleep_sub_stage_cnt, 2)
        if sleep_sub_stage_cnt(stagei, subi) == 0
            continue;
        end
        nc_time = nc_time+1;
        Modularity_all_cell(nc_time, :) = [M_all_list(subi, stagei) stagei AgeInYears(subi) Gender(subi) EducationYear(subi) subi sleep_sub_stage_cnt(stagei, subi) FD_stage_list(stagei, subi)];
    end
end
t1 = array2table(Modularity_all_cell, 'VariableNames', variable_Names_list);
writetable(t1, [ ret_dir '/Modularity_all_table_individual_atlas_cluster.csv']);
