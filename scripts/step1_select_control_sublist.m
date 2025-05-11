% Author: Jinlong Li and Guoyuan Yang, BIT.
% select control group for validation
% 

clc;clear;

work_dir = '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod';
ret_dir = '/nd_disk3/guoyuan/Jinlong/sleep_ret/HFR_ret';


work_dir = work_code_dir;
ret_dir = sleep_ret_dir;

% network_label_all sleep_sub_file_info sleep_sub_stage_cnt sublist sleep_sub_stage_info  filter_sub_idx 
load([ ret_dir '/HFR_ret/IndiPar_net_multi_7/Network_ind_subcortex_cleanup5min.mat']);
load([ work_dir '/resource/fs_LR_32k_no_medial_wall_index.mat']);  %no_medial_wall_index




load([ work_dir '/subject_info.mat']);

min_subn = min(sum(sleep_sub_stage_cnt(1:4, :) ~= 0, 2));
min_subn = floor(min_subn / 2)*2;


network_n = 7;
vertex_n = 96854;



% select subjects set with balanced FD (Mannual, hard to select randomly for FD)
p_list_cnt_total = cell(22, 1);
control_sub_list = {}; 
for rand_i = 1:50
    stage_select_list = { ...
        [53-rand_i:52 53+rand_i:128],...
        [1:76],...
        [1:76],...
        [1:76]}; 
    p_list_tmp = zeros(4, 4);
    p_flg = 0;
    control_sub_idx = cell(4, 1);
    for stagei = 1:3
        for stagej = (stagei+1) : 4
            ka = FD_stage_list(stagei, :);
            ka(ka == 0) = [];
            [kv1, ki1] = sort(ka);
            ki2 = ki1(stage_select_list{stagei});
            kc = ka(ki2);
            tmp1 = find(FD_stage_list(stagei, :) ~= 0);
            control_sub_idx{stagei} = tmp1(ki2);

            ja = FD_stage_list(stagej, :);
            ja(ja == 0) = [];
            [jv1, ji1] = sort(ja);
            ji2 = ji1(stage_select_list{stagej});
            jc = ja(ji2);
            tmp1 = find(FD_stage_list(stagej, :) ~= 0);
            control_sub_idx{stagej} = tmp1(ji2);

            [h, p] = ttest2(kc, jc);
            p_list_tmp(stagei, stagej) = p;
            if p <= 0.05
                p_flg = 1;
            end
        end
    end
    p_list_cnt_total{rand_i} = p_list_tmp;
    if p_flg == 0
        control_sub_list{end+1} = control_sub_idx;
    end
end



control_sub_list_select = control_sub_list{1};

save([work_dir '/subject_control_select_7.mat'], 'control_sub_list', 'control_sub_list_select', 'p_list_cnt_total');