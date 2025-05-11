% Author: Jinlong Li and Guoyuan Yang, BIT.
% cleanup for subcortex atlas constructed with all sessions
% 

clc;clear;

work_dir = '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/';
ret_dir = '/nd_disk3/guoyuan/Jinlong/sleep_ret/HFR_ret/';

addpath '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/code'

% network_label_all sleep_info sublist sub_stage_cell 
load([ret_dir '/IndiPar_net7/Network_ind_ret8min.mat']);
load([ work_dir '/resource/fs_LR_32k_no_medial_wall_index.mat']);  %no_medial_wall_index
load([ work_dir '/sleep_subinfo.mat']);
% line 6: eo1
sleep_sub_stage_cnt(:, ~filter_sub_idx) = [];
sleep_sub_file_info(:, ~filter_sub_idx) = [];
sleep_sub_stage_info(:, ~filter_sub_idx) = [];
sublist(~filter_sub_idx) = [];
network_label_all(:, ~filter_sub_idx) = [];

stageNames = {'Awake', 'N1', 'N2', 'N3', 'REM'};

load([ret_dir '/IndiPar_net7/Network_ind_subcortex8min.mat']);


%% clean up
template_dlabel = '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/resource/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_netassignments_LR.dlabel.nii';
template_MNI152_cortex = '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/resource/MNI152_cortex_mask.nii.gz';

%% dist_subcortex_subcortex
dlabel_struct = ft_read_cifti(template_dlabel);
dist_subcortex_subcortex = compute_vertices_dist(dlabel_struct.pos(64985:end, :), dlabel_struct.pos(64985:end, :));

%% neighbor_subcortex
neighbor_subcortex = cell(96854-64984, 1);
for voxel_i = 1:(96854-64984)
    neighbor_list = find_neighbor_voxel(voxel_i, dist_subcortex_subcortex, dlabel_struct.pos(64985:end, :));
    neighbor_subcortex{voxel_i} = neighbor_list;
end

%% cerebellum_voxel_idx
MNI_ras2vox_file(1:96854, './test.nii.gz');
nifti_base_info = MRIread('./test.nii.gz');
system('rm ./test.nii.gz');
template_MNI152_cortex = '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/resource/MNI152_cortex_mask.nii.gz';
MNI152_struct = MRIread(template_MNI152_cortex);
nzk = MNI152_struct.vol ~= 0 & nifti_base_info.vol ~= 0;
overlap_idx = nifti_base_info.vol(nzk);
cerebellum_voxel_idx = zeros(96854, 1);
cerebellum_voxel_idx(overlap_idx) = 1;



%% cleanup

subcortex_label = network_label_subcortex_lisa;
for stagei = 1:5
    for si = 1:size(subcortex_label, 2)
        disp(['process stage:' num2str(stagei) ' subi:' num2str(si) ]);tic;
        if isempty(subcortex_label{stagei, si})
            continue;
        end
        label1 = subcortex_label{stagei, si};
        new_label_assign = subcortex_parcellations_cleanup(label1, 'cerebellum_voxel_idx', cerebellum_voxel_idx, 'dist_subcortex_subcortex', dist_subcortex_subcortex, 'neighbor_subcortex', neighbor_subcortex);
        disp(['change ratio:' num2str(nnz(new_label_assign ~= subcortex_label{stagei, si}) / 96854)]);
        subcortex_label{stagei, si} = new_label_assign;
        
        toc;
    end
end
network_label_subcortex_lisa = subcortex_label;

