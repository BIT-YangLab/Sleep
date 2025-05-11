% Author: Jinlong Li and Guoyuan Yang, BIT.
% prepare information for analysis 
% 


clc;clear;

work_dir = '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod';
ret_dir = '/nd_disk3/guoyuan/Jinlong/sleep_ret/HFR_ret';



Ancillary_data_file = [ work_dir '/resource/Ancillary_data_final.xlsx'];
Model_sleep_data_file = [work_dir '/resource/Model_sleep_wake_SWA_N2562.txt'];


% network_label_all sleep_info sublist sub_stage_cell 
% line 6: EO as awake stage
load([ work_dir '/sleep_subinfo.mat']);
sleep_sub_stage_cnt(:, ~filter_sub_idx) = []; sleep_sub_stage_cnt = sleep_sub_stage_cnt([6 2 3 4 5], :);
sleep_sub_file_info(:, ~filter_sub_idx) = []; sleep_sub_file_info = sleep_sub_file_info([6 2 3 4 5], :);
sleep_sub_stage_info(:, ~filter_sub_idx) = []; sleep_sub_stage_info = sleep_sub_stage_info([6 2 3 4 5], :);
sublist(~filter_sub_idx) = []; 


load([ work_dir '/resource/fs_LR_32k_no_medial_wall_index.mat']);  %no_medial_wall_index

stageNames = {'Awake', 'N1', 'N2', 'N3', 'REM'};



sub_number = sum(sleep_sub_stage_cnt ~= 0, 2);
ses_number = sum(sleep_sub_stage_cnt, 2);

Ancillary_data = readtable(Ancillary_data_file);
Model_sleep_data = readtable(Model_sleep_data_file);

subID = Ancillary_data.ID;

sublist1 = zeros(size(subID, 2), 1);
for ki = 1:size(subID, 1)
    tmp_str = strsplit(subID{ki}, 'sub');
    sublist1(ki) = str2double(tmp_str{2});
end

[ia, ib] = ismember(sublist, sublist1);

Ancillary_data = Ancillary_data(ib, :);

% age distribution
age = Ancillary_data.AgeInYears;
age_list = cell(5, 1);
for stagei = 1:5
    idx = find(sleep_sub_stage_cnt(stagei, :) ~= 0);
    age_list{stagei} = age(idx);
    disp(['stage: ' stageNames{stagei} ' age range: ' num2str(min(age_list{stagei})) '~' num2str(max(age_list{stagei})) ... 
     ' mean:' num2str(mean(age_list{stagei})) ' std:' num2str(std(age_list{stagei}))]);
end

% sex distribution 0-female 1-male 
sex = Ancillary_data.Gender;
sex_list = zeros(5, 2);
for stagei = 1:5
    idx = find(sleep_sub_stage_cnt(stagei, :) ~= 0);
    sex_list(stagei, 1) = nnz(sex(idx) == 0) ./ length(idx);
    sex_list(stagei, 2) = nnz(sex(idx)) ./ length(idx);
    disp(['stage:' stageNames{stagei} ' sex: female-' num2str(sex_list(stagei, 1)) ' male-' num2str(sex_list(stagei, 2)) ]);
end

DV_cell = cell(stagei, 1); FD_cell = cell(stagei, 1);
DV_stage_cell = cell(size(sleep_sub_stage_cnt));
FD_stage_cell = cell(size(sleep_sub_stage_cnt));
DV_stage_list = zeros(size(sleep_sub_stage_cnt));
FD_stage_list = zeros(size(sleep_sub_stage_cnt));
for stagei = 1:5
    DV_list = []; FD_list = [];
    for subi = 1:size(sleep_sub_stage_info, 2)
        if isempty(sleep_sub_stage_info{stagei, subi})
            continue;
        end
        sesn = sleep_sub_stage_info{stagei, subi};
        subname = num2str(sublist(subi));
        DV_list_tmp = []; FD_list_tmp = []; DV_list_mean_tmp = []; FD_list_mean_tmp = [];
        for sesi = 1:length(sesn)
            if stagei == 1
                sess = '01';
            else
                sess = sesn{sesi};
            end
            if stagei == 1 && sesi == 1
                DV_FD_file = ['/nd_disk3/guoyuan/sleep/sleep_eo_ciftify_1/out/fmriprep/sub-' subname '/ses-' sess '/func/sub-' subname '_ses-' sess '_task-sleep_desc-confounds_regressors.tsv'];
            elseif stagei == 1 && sesi == 2
                DV_FD_file = ['/nd_disk3/guoyuan/sleep/sleep_eo_ciftify_2/out/fmriprep/sub-' subname '/ses-' sess '/func/sub-' subname '_ses-' sess '_task-sleep_desc-confounds_regressors.tsv'];
            else 
                DV_FD_file = ['/md_disk4/guoyuan/Sleep/Sleep_preproc/Sleep_out_preproc_300s/out/fmriprep/sub-' subname '/ses-' sess '/func/sub-' subname '_ses-' sess '_task-sleep_desc-confounds_regressors.tsv'];
            end
            DV_FD_info = readtable(DV_FD_file, 'FileType', 'text');
            FD = DV_FD_info.framewise_displacement;
            DV = DV_FD_info.dvars;
            FD_list = [ FD_list; FD]; DV_list = [ DV_list; DV];
            FD_list_mean_tmp = [FD_list_mean_tmp; nanmean(FD)]; DV_list_mean_tmp = [DV_list_mean_tmp; nanmean(DV)];
        end
        FD_stage_cell{stagei, subi} = FD_list_mean_tmp; FD_stage_list(stagei, subi) = nanmean(FD_list);
        DV_stage_cell{stagei, subi} = DV_list_mean_tmp; DV_stage_list(stagei, subi) = nanmean(DV_list);
    end
    FD_list(isnan(FD_list)) = [];
    DV_list(isnan(DV_list)) = [];
    FD_cell{stagei} = FD_list;
    DV_cell{stagei} = DV_list;
end

for stagei = 1:5
    disp(['stage:' stageNames{stagei} ' FD: mean-' num2str(nanmean(FD_cell{stagei})) ' std-' num2str(std(FD_cell{stagei}))]);
end

for stagei = 1:5
    disp(['stage:' stageNames{stagei} ' DV: mean-' num2str(nanmean(DV_cell{stagei})) ' std-' num2str(std(DV_cell{stagei}))]);
end

PSQI = Ancillary_data.PSQI;
ESS = Ancillary_data.ESS;
MEQ = Ancillary_data.MEQ;
FSS = Ancillary_data.FSS;

% SWA
SWA = Model_sleep_data.SWA;
Model_file = Model_sleep_data.filename;
swa_cell = cell(size(sleep_sub_stage_info));

format_EEG_file = '/md_disk4/guoyuan/Sleep/Sleep_preproc/Sleep_out_preproc_300s/out/filename_txt/filename_%d.txt';


for stagei = 2:5
    for subi = 1:size(sleep_sub_stage_info, 2)
        swa_list_tmp = [];
        if isempty(sleep_sub_stage_info{stagei, subi})
            swa_cell{stagei, subi} = swa_list_tmp;
            continue;
        end
        sesn = sleep_sub_stage_info{stagei, subi};
        EEG_stage_name = importdata(sprintf(format_EEG_file, sublist(subi)));
        
        for sesi = 1:length(sesn)
            ses_tmp = sesn{sesi};
            str_idx = sprintf('sub%d_%s_1', sublist(subi), EEG_stage_name{str2double(ses_tmp)});
            kidx = find(strcmp(Model_file, str_idx));
            if isempty(kidx)
                swa_list_tmp(end+1) = nan;
            else
                swa_list_tmp(end+1) = SWA(kidx);
            end
            
        end
        swa_cell{stagei, subi} = swa_list_tmp;
    end
end

el_value = Model_sleep_data.sw_history_2bins;
early_late_cell = cell(size(sleep_sub_stage_info));
for stagei = 2:5
    for subi = 1:size(sleep_sub_stage_info, 2)
        el_list_tmp = [];
        if isempty(sleep_sub_stage_info{stagei, subi})
            early_late_cell{stagei, subi} = el_list_tmp;
            continue;
        end
        sesn = sleep_sub_stage_info{stagei, subi};
        EEG_stage_name = importdata(sprintf(format_EEG_file, sublist(subi)));
        
        for sesi = 1:length(sesn)
            ses_tmp = sesn{sesi};
            str_idx = sprintf('sub%d_%s_1', sublist(subi), EEG_stage_name{str2double(ses_tmp)});
            kidx = find(strcmp(Model_file, str_idx));
            if isempty(kidx)
                el_list_tmp(end+1) = nan;
            else
                el_list_tmp(end+1) = strcmp(el_value{kidx}, 'onemore');
            end
            
        end
        early_late_cell{stagei, subi} = el_list_tmp;
    end
end

save([ work_code_dir '/subject_info.mat' ], ...
    'sleep_sub_stage_cnt', 'sleep_sub_file_info', 'sleep_sub_stage_info', ...
    'sublist', 'Ancillary_data', 'Model_sleep_data', ...
    'age_list', 'sex_list', 'DV_cell', 'FD_cell', 'swa_cell', ...
    'DV_stage_cell', 'FD_stage_cell', 'DV_stage_list', 'FD_stage_list');




