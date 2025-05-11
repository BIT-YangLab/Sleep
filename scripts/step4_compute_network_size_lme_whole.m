% Author: Jinlong Li and Guoyuan Yang, BIT.
% compute individual network size 
% 
clc;clear;

work_dir = '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/';
ret_dir = '/nd_disk3/guoyuan/Jinlong/sleep_ret/HFR_ret/';

% network_label_all sleep_info sublist sub_stage_cell 
load([ret_dir '/IndiPar_net_multi_7/Network_ind_ret5min.mat']);
load('/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/sleep_subinfo.mat');
load('/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/resource/fs_LR_32k_no_medial_wall_index.mat');  %no_medial_wall_index
% line 6: eo1
sleep_sub_stage_cnt(:, ~filter_sub_idx) = [];
sleep_sub_file_info(:, ~filter_sub_idx) = [];
sleep_sub_stage_info(:, ~filter_sub_idx) = [];
sublist(~filter_sub_idx) = [];
network_label_all(:, ~filter_sub_idx) = [];

network_n = 7;
vertex_n = 64984;

variable_Names_list = {'NS', 'Stage', 'Age', 'Gender', 'Education', 'participant_id', 'epoch', 'fd'};

network_size_ret = cell(size(sleep_sub_stage_info));


load([ret_dir '/IndiPar_net7/Network_ind_subcortex_cleanup8min.mat']);

template_dlabel = '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/resource/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_netassignments_LR.dlabel.nii';

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
mask_list_w = {cortex_mask_lh | cortex_mask_rh, cerebellum_mask_lh | cerebellum_mask_rh, thalamus_mask_lh | thalamus_mask_rh, striatum_mask_lh | striatum_mask_rh};
mask_name_list_w = {'cortex', 'cerebellum', 'thalamus', 'striatum'};

group_label_list = cell(4, 1);

load([work_dir '/subject_info.mat']);

AgeInYears = Ancillary_data.AgeInYears;
Gender = Ancillary_data.Gender;
EducationYear = Ancillary_data.EducationInYears;

subcortex_label = network_label_subcortex_lisa;
label1 = [];


network_label_all = network_label_subcortex_lisa;

network_size_all_list = cell(length(mask_list_w), 1);

for mask_i = 1:length(mask_list_w)
    for stai = 1:5
        
        nc = 0;
        for subi = 1:length(network_label_all)
            if isempty(network_label_all{stai, subi})
                continue;
            end
    
            temp_size_sta = zeros(1, network_n);
            
            network_tmp_label = network_label_all{stai, subi};
            network_tmp_label = network_tmp_label(mask_list_w{mask_i});
            for ni = 1:network_n
                idx = find(network_tmp_label == ni);
                if isempty(idx)
                    disp([num2str(stai) ' ' num2str(subi)])
                end
                temp_size_sta(1, ni) = length(idx);
            end
           
            network_size_ret{stai, subi} = temp_size_sta./length(network_tmp_label);
        end
        
    end

    network_size_all_cell = cell(network_n, 1);
    nc_time = 0;
    for stagei = 1:5
        for subi = 1:length(sublist)
            if isempty(network_size_ret{stagei, subi}) 
                continue;
            end
            nc_time = nc_time + 1;
            for networki = 1:network_n
                if isempty(network_size_all_cell{networki}) 
                    network_size_all_cell{networki} = zeros(10, 8);
                end
                network_size_all_cell{networki}(nc_time, :) = [network_size_ret{stagei, subi}(networki) stagei AgeInYears(subi) Gender(subi) EducationYear(subi) subi sleep_sub_stage_cnt(stagei, subi) FD_stage_list(stagei, subi)];
            end
        end
    end


    network_size_cell = cell(network_n, 1);
    for stai = 1:5
        nc = 0;
        for subi = 1:length(sublist)
            if isempty(network_size_ret{stai, subi})
                continue;
            end 
            nc = nc + 1;
            for neti = 1:network_n
                if isempty(network_size_cell{neti})
                    network_size_cell{neti} = zeros(10, 4);
                end
                
                network_size_cell{neti}(nc, stai) = network_size_ret{stai, subi}(neti);
            end

        end
    end
    
    network_size_cell2 = cell(network_n, 1);


    p_list_cell = cell(network_n, 1);
    
    stageNames = {'Awake', 'N1', 'N2', 'N3', 'REM'};
    
    % write table into file and compute the significance by LME in Rscript
    for neti = 1:network_n
        temp_net_size = cell(6, 1);
        temp_ns = network_size_cell{neti};
        temp_ns(temp_ns == 0) = nan;
        network_size_cell{neti} = temp_ns;

        t1 = array2table(network_size_all_cell{neti, 1}, 'VariableNames', variable_Names_list);
        writetable(t1, ['/nd_disk3/guoyuan/Jinlong/sleep_ret/HFR_ret/Network_size_table_brainstructure_' num2str(mask_i) '_net' num2str(neti) '.csv']);
    end


    network_size_all_list{mask_i} = network_size_cell;


end







