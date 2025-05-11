% Author: Jinlong Li and Guoyuan Yang, BIT.
% construct density map and group average atlas 
% 

clc;clear;

work_dir = '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai1_mod/';
ret_dir = '/nd_disk3/guoyuan/Jinlong/sleep_ret/HFR_ret/';

% network_label_all sleep_info sublist sub_stage_cell 
load([ret_dir '/IndiPar_net7/Network_ind_ret8min.mat']);
load([work_dir '/sleep_subinfo.mat']);
load([work_dir '/resource/fs_LR_32k_no_medial_wall_index.mat']);  %no_medial_wall_index


addpath(['work_dir '/code']);



% line 6: eo1
sleep_sub_stage_cnt(:, ~filter_sub_idx) = [];
sleep_sub_file_info(:, ~filter_sub_idx) = [];
sleep_sub_stage_info(:, ~filter_sub_idx) = [];
sublist(~filter_sub_idx) = [];
network_label_all(:, ~filter_sub_idx) = [];

network_n = 7;
vertex_n = 64984;

network_ret = cell(6, 1);

for stai = 1:6
    temp_label_sta = zeros(vertex_n, network_n);
    nc = 0;
    for subi = 1:length(sublist)
        if isempty(network_label_all{stai, subi})
            continue;
        end
        for ni = 1:network_n
            idx = find(network_label_all{stai, subi} == ni);
            temp_label_sta(idx, ni) = temp_label_sta(idx, ni) + 1;
        end
        nc = nc + 1;
    end
    temp_label_sta = temp_label_sta ./ nc;
    network_ret{stai} = temp_label_sta;
end

dt = ft_read_cifti('/nd_disk3/guoyuan/Jinlong/test.dtseries.nii');
dt.time = 1:network_n;

%% network
group_density_map = zeros(vertex_n, 6);
for stai = 1:6
    temp_label_sta = network_ret{stai};
    temp_label_sta(temp_label_sta < 0.4) = 0;
    temp_label_sta(temp_label_sta > 0.8) = 0.8;
    dt.dtseries = temp_label_sta;
    group_density_map(:, stai) = temp_label_sta(:, 4);
    ft_write_cifti([ret_dir '/sleep_st' num2str(stai)  ], dt, 'parameter', 'dtseries');
end

%% whole brain (winner-take-all)
group_label = zeros(vertex_n, 6);
for stai = 1:6
    temp_label_sta = network_ret{stai};
    [~, colIndex] = max(temp_label_sta, [], 2);
    colIndex(setdiff(1:vertex_n, no_medial_wall_index)) = 0;
    group_label(:, stai) = colIndex;

end


 


