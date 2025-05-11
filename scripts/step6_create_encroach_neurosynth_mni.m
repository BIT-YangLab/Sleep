% Author: Jinlong Li and Guoyuan Yang, BIT.
% select reduction region and transfer data into MNI152 2mm space
% get functional terms from the results from Neurosynth
clc;clear;

setenv('CBIG_CODE_DIR', '/nd_disk3/guoyuan/Jinlong/CBIG');
addpath('/nd_disk3/guoyuan/Jinlong/CBIG/utilities/matlab/transforms');


stage_names = {'Awake', 'N1', 'N2', 'N3', 'N4'};

% neurosynth
load('/nd_disk3/guoyuan/Jinlong/sleep_ret/HFR_ret/IndiPar_net7/Network_grp_amn_density_ret.mat')
load('/nd_disk3/guoyuan/Jinlong/sleep_ret/HFR_ret/IndiPar_net7/Network_grp_ret.mat')
g_amn = group_label(:, 1);
g_amn(g_amn ~= 4) = 0;

for stagei = 1:4
    a1 = group_label(:, stagei);
    a1(a1 ~= 4) = 0;
    color_label_visual(a1, '/nd_disk3/guoyuan/Jinlong/sleep_ret/HFR_ret/', ['_encroach_' stage_names{stagei} '_'], 64984, 7);

end

% sleep stage
for stagei = 2:4
    g_n1 = group_label(:, stagei);
    g_n1(g_n1 ~= 4) = 0;
    idx1 = find(g_amn ~= 0);
    idx2 = find(g_n1 ~= 0);
    idx3 = setdiff(idx1, idx2);
    grp_df = group_label(:, stagei);
    grp_df(idx3) = 8;
    [lh_FS7_data,rh_FS7_data] = CBIG_project_fsLR2fsaverage(grp_df(1:32492, 1), grp_df(32493:64984, 1), 'fs_LR_32k', 'label', '/nd_disk3/guoyuan/Jinlong/sleep_ret/HFR_ret/sleep_density_temp');
    output = CBIG_Projectfsaverage2MNI_Ants(lh_FS7_data, rh_FS7_data, 'nearest');
    MRIwrite(output, '/nd_disk3/guoyuan/Jinlong/sleep_ret/MNI_256_label.nii.gz');
    system(['flirt -in /nd_disk3/guoyuan/Jinlong/sleep_ret/MNI_256_label.nii.gz -ref /nd_disk3/guoyuan/Jinlong/sleep_ret/MNI152_T1_2mm_brain.nii.gz -out /nd_disk3/guoyuan/Jinlong/sleep_ret/HFR_ret/MNI_stage_' stage_names{stagei} '_encroach.nii.gz -interp nearestneighbour;'])
end

%% after computing neurosynth 
n1 = readtable('/nd_disk3/guoyuan/Jinlong/sleep_ret/HFR_ret/decoding_results_N1.txt');
n2 = readtable('/nd_disk3/guoyuan/Jinlong/sleep_ret/HFR_ret/decoding_results_N2.txt');
n3 = readtable('/nd_disk3/guoyuan/Jinlong/sleep_ret/HFR_ret/decoding_results_N3.txt');
[iv, ia] = sort(n1.cluster_7, 'descend');
feature_n1 = n1.Feature(ia);
meta_N1 = table(feature_n1, iv, 'VariableNames', {'Feature', 'zscore'});
[iv, ia] = sort(n2.cluster_7, 'descend');
feature_n2 = n2.Feature(ia);
meta_N2 = table(feature_n2, iv, 'VariableNames', {'Feature', 'zscore'});
[iv, ia] = sort(n3.cluster_7, 'descend');
feature_n3 = n3.Feature(ia);
meta_N3 = table(feature_n3, iv, 'VariableNames', {'Feature', 'zscore'});

% only care about Zscore > 0
Z = meta_N1.zscore; Z(Z <= 0) = [];
p_value = 1 - normcdf(Z);
p_value_l = meta_N1.zscore(meta_N1.zscore <= 0);
p_value_l(:) = 1;
p_value = [p_value; p_value_l];
[pthr,pcor,q_list] = fdr(p_value);
p_list_N1 = q_list;

Z = meta_N2.zscore; Z(Z <= 0) = [];
p_value = 1 - normcdf(Z);
p_value_l = meta_N2.zscore(meta_N2.zscore <= 0);
p_value_l(:) = 1;
p_value = [p_value; p_value_l];
[pthr,pcor,q_list] = fdr(p_value);
p_list_N2 = q_list;

Z = meta_N3.zscore; Z(Z <= 0) = [];
p_value = 1 - normcdf(Z);
p_value_l = meta_N3.zscore(meta_N3.zscore <= 0);
p_value_l(:) = 1;
p_value = [p_value; p_value_l];
[pthr,pcor,q_list] = fdr(p_value);
p_list_N3 = q_list;
