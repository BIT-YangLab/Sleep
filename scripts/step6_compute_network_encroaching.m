
clc;clear;

work_dir = '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/';
ret_dir = '/nd_disk3/guoyuan/Jinlong/sleep_ret/HFR_ret/';

% network_label_all sleep_info sublist sub_stage_cell 
load([ret_dir '/IndiPar_net7/Network_ind_ret8min.mat']);
% group_label
load([ret_dir '/IndiPar_net7/Network_grp_ret.mat']);
load([ work_dir '/sleep_subinfo.mat']);
load([ work_dir '/resource/fs_LR_32k_no_medial_wall_index.mat']);  %no_medial_wall_index
% line 6: eo1
sleep_sub_stage_cnt(:, ~filter_sub_idx) = [];
sleep_sub_file_info(:, ~filter_sub_idx) = [];
sleep_sub_stage_info(:, ~filter_sub_idx) = [];
sublist(~filter_sub_idx) = [];
network_label_all(:, ~filter_sub_idx) = [];

network_n = 7;
stageNames = {'Awake', 'N1', 'N2', 'N3', 'REM'};
% ventral attention
network_specific = 4;

relist = [6 2 3 4 5];
network_label_all = network_label_all(relist, :);
% group_label = group_label(:, relist);

% encroaching prepare
dt = ft_read_cifti('/nd_disk3/guoyuan/Jinlong/test.dtseries.nii');
dt.time = 1;
mid_L = '/nd_disk3/guoyuan/Jinlong/Gordon2016_2/Parcels/Conte69.L.midthickness.32k_fs_LR.surf.gii';
mid_R = '/nd_disk3/guoyuan/Jinlong/Gordon2016_2/Parcels/Conte69.R.midthickness.32k_fs_LR.surf.gii';
encraoaching_cluster_file_stem = '/nd_disk3/guoyuan/Jinlong/sleep_ret/HFR_ret/Temp_encroaching_cluster';
sizeThresh = 30;
distThresh = 3.5;
load('/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/resource/fs_LR_32k_surf_distance_lh.mat');
load('/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/resource/fs_LR_32k_surf_distance_rh.mat');
dist_matrix = zeros(64984, 64984);
dist_matrix(1:32492, 1:32492) = lh_dist;
dist_matrix(32493:64984, 32493:64984) = rh_dist;
dist_matrix(1:32492, 32493:64984) = 255;
dist_matrix(32493:64984, 1:32492) = 255;
clear lh_dist rh_dist

encroaching_ret = cell(size(network_label_all));
encroaching_type_list = cell(5, 1);
encroaching_network_list = cell(5, 1);

variable_Names_list = {'NS', 'Stage', 'Age', 'Gender', 'Education', 'participant_id', 'epoch', 'fd'};

load([work_dir '/subject_info.mat']);

AgeInYears = Ancillary_data.AgeInYears;
Gender = Ancillary_data.Gender;
EducationYear = Ancillary_data.EducationInYears;


network_size_all_cell = cell(4, 1);
for stagei = 2:4
    encroaching_type_tmp = zeros(10, 2, 2);
    encroaching_network_tmp = zeros(10, 7, 2);
    nc = 0;
    nc_time = 0;
    for subi = 1:size(network_label_all, 2)
        if isempty(network_label_all{stagei, subi})
            continue;
        end
        disp(['process stage:' num2str(stagei) ' subi:' num2str(subi) '/' num2str(size(network_label_all, 2)) ]); tic;

%         label_1 = group_label(:, 1); label_1(label_1 ~= network_specific) = 0; 
        label_1 = network_label_all{1, subi}; label_1(label_1 ~= network_specific) = 0;
        label_2 = network_label_all{stagei, subi}; label_2(label_2 ~= network_specific) = 0;

        % cluster group_label(mask) encroach_Type(0-intrusion)
        variant_cluster_label = zeros(size(label_1, 1), 3, 2);

        % q. decrease or increase
        label_diff = (label_1 - label_2); label_diff(64985:96854) = 0;
        label_intersect = (label_1 + label_2) == 2*network_specific;
        nc = nc + 1;
        for ti = 1:2
            if ti == 1
                dt.dtseries = label_diff > 0;
                template_label = network_label_all{stagei, subi};
            else
                dt.dtseries = label_diff < 0;
                template_label = group_label(:, 1);
            end
            ft_write_cifti(encraoaching_cluster_file_stem, dt, 'parameter', 'dtseries')
            system(['wb_command -cifti-find-clusters ' encraoaching_cluster_file_stem '.dtseries.nii  0 0 0 0 COLUMN ' encraoaching_cluster_file_stem '_ret.dtseries.nii -left-surface ' ... 
            mid_L ' -right-surface ' mid_R ]);
            cluster_dt = ft_read_cifti([encraoaching_cluster_file_stem '_ret.dtseries.nii']);
            new_label = cluster_dt.dtseries(1:64984);
            
            cluster_l = unique(new_label); cluster_l(cluster_l < 1) = [];
            
            for cli = 1:length(cluster_l)
                if nnz(new_label == cluster_l(cli)) < sizeThresh
                    new_label(new_label == cluster_l(cli)) = 0;
                end
            end
            
            variant_cluster_label(:, 1, ti) = new_label;
            variant_cluster_label(new_label ~= 0, 2, ti) = template_label(new_label ~= 0);

            % ectopic and intrusion
            cluster_l = unique(new_label); cluster_l(cluster_l < 1) = [];
            cluster_number_list = zeros(length(cluster_l), 1);
            group_network_idx = group_label(:, 1) == network_specific;
            for cli = 1:length(cluster_l)
                idx_tmp = new_label == cluster_l(cli);
                dist_tmp = dist_matrix(idx_tmp, label_intersect);
                min_dist = min(min(dist_tmp));
                variant_cluster_label(idx_tmp, 3, ti) = min_dist > distThresh;
            end
            
            encroaching_type_tmp(nc, 1, ti) = nnz((variant_cluster_label(:, 1, ti) ~= 0) - variant_cluster_label(:, 3, ti)) / nnz(variant_cluster_label(:, 1, ti));
            encroaching_type_tmp(nc, 2, ti) = nnz(variant_cluster_label(:, 3, ti)) / nnz(variant_cluster_label(:, 1, ti));
            
            for networki = 1:network_n
                encroaching_network_tmp(nc, networki, ti) = nnz(variant_cluster_label(:, 2, ti) == networki) / nnz(variant_cluster_label(:, 1, ti));
            end
            
        end
        encroaching_ret{stagei, subi} = variant_cluster_label;

        if isempty(network_size_all_cell{stagei}) 
            network_size_all_cell{stagei} = zeros(10, 8);
        end
        for networki = 1:network_n
            if networki == network_specific
                continue;
            end
            nc_time = nc_time + 1;
            network_size_all_cell{stagei}(nc_time, :) = [encroaching_network_tmp(nc, networki, 1) networki AgeInYears(subi) Gender(subi) EducationYear(subi) subi sleep_sub_stage_cnt(stagei, subi) FD_stage_list(stagei, subi)];
        end

        toc;
    end
    encroaching_type_list{stagei} = encroaching_type_tmp;   
    encroaching_network_list{stagei} = encroaching_network_tmp;

    t1 = array2table(network_size_all_cell{stagei, 1}, 'VariableNames', variable_Names_list);
    writetable(t1, ['/nd_disk3/guoyuan/Jinlong/sleep_ret/HFR_ret/text/Encoaching_table_stage' stageNames{stagei} '.csv']);

    
end


% intrusion(idx=0) and ectopic(idx=1) - encroaching type
% reduction(ti=1) and expansion(ti=2) 


ecr_type_list = cell(5, 1);
for stagei = 1:5
    nc = 0;
    ecr_type = zeros(10, 2);
    for subi = 1:size(encroaching_ret, 2)
        if isempty(encroaching_ret{stagei, subi})
            continue;
        end
        nc = nc + 1;
        ecr_type(nc, 1) = nnz(encroaching_ret{stagei, subi}(:, 1, 1)) ./ 59412;
        ecr_type(nc, 2) = nnz(encroaching_ret{stagei, subi}(:, 1, 2)) ./ 59412;
    end
    ecr_type_list{stagei} = ecr_type;
end
