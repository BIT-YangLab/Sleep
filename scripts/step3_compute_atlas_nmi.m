% Author: Jinlong Li and Guoyuan Yang, BIT.
% compute nmi of group atlas between control group and primary analysis group 
% 

clc;clear;

work_dir = '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/';
ret_dir = '/nd_disk3/guoyuan/Jinlong/sleep_ret/HFR_ret/';

addpath '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/code'

% network_label_all sleep_info sublist sub_stage_cell 

load([ work_dir '/resource/fs_LR_32k_no_medial_wall_index.mat']);  %no_medial_wall_index
load([ work_dir '/sleep_subinfo.mat']);



stageNames = {'Awake', 'N1', 'N2', 'N3', 'REM'};

load([ret_dir '/IndiPar_net7/Network_ind_subcortex_cleanup8min.mat']);



%% clean up
template_dlabel = '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/resource/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_netassignments_LR.dlabel.nii';
template_MNI152_cortex = '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/resource/MNI152_cortex_mask.nii.gz';


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
mask_list_w = {cortex_mask_lh | cortex_mask_rh, cerebellum_mask_lh | cerebellum_mask_rh, thalamus_mask_lh | thalamus_mask_rh, striatum_mask_lh | striatum_mask_rh};
mask_name_list_w = {'cortex', 'cerebellum', 'thalamus', 'striatum'};


load([ret_dir '/IndiPar_net7/group_atlas.mat']);
ga = group_label_list;
load([ret_dir '/IndiPar_net_multi_7/group_control_atlas.mat']);
gc = group_label_list;




nmi_list_control = zeros(4, 4);


dice_list_control = zeros(4, 4);


for stagei = 1:4
    l1 = ga{stagei}; l2 = gc{stagei};
    for maski = 1:4
        l_tmp1 = l1(mask_list_w{maski}); l_tmp2 = l2(mask_list_w{maski});
        nmi_value = nmi(l_tmp1, l_tmp2);
        nmi_list_control(stagei, maski) = nmi_value;

        [dice_list, ~] = compute_dice_coef(l_tmp1, l_tmp2);
        
        dice_list_control(stagei, maski) = mean(dice_list);
    end
end
