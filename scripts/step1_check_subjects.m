% Author: Jinlong Li and Guoyuan Yang, BIT.
% check subject info from sleep datasets
% 


clc;clear;

work_code_dir = '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod';
sleep_ret_dir = '/nd_disk3/guoyuan/Jinlong/sleep_ret/HFR_ret';

Sleep_300s_post_file = '/nd_disk3/guoyuan/sleep/sleep_300s_ciftify_post_new/postpreprocess_ciftify/';
Sleep_300s_pre_file = '/md_disk4/guoyuan/Sleep/Sleep_preproc/Sleep_out_preproc_300s/out/'; 
EO_post_file = '/nd_disk3/guoyuan/sleep/sleep_eo_ciftify_%d/postpreprocess_ciftify/';
Sleep_300s_regress_file = '/md_disk4/guoyuan/Sleep/Sleep_preproc/Sleep_out_preproc_300s/out/fmriprep/sub-%d/ses-%s/func/sub-%d_ses-%s_task-sleep_desc-confounds_regressors.tsv';
EO_post_regress_file = '/nd_disk3/guoyuan/sleep/sleep_eo_ciftify_%s/out/fmriprep/sub-%d/ses-0%s/func/sub-%d_ses-0%s_task-sleep_desc-confounds_regressors.tsv';


sub_list_file = [ work_code_dir 'sleep_subjects_list.txt'];

sublist = importdata(sub_list_file);

% awake, N1, N2, N3, REM, EO (eyes open)
sleep_sub_stage_info = cell(6, length(sublist));
sleep_sub_file_info = cell(6, length(sublist));
sleep_sub_stage_cnt = zeros(6, length(sublist));



for subi = 1:length(sublist)
    subname = num2str(sublist(subi));
    file1 = [ Sleep_300s_pre_file '/filename_txt/filename_' subname '.txt'];
    if ~exist(file1)
        continue;
    end
    txt_line = importdata(file1);
    for li = 1:length(txt_line)
        line1 = txt_line{li}; % stage{stagei}sess{sessi}
        stagei = sscanf(line1,'stage%dsess%d');
        stagei = stagei(1);
        sess_format = sprintf('%02d', li);
        cifti_file = [ Sleep_300s_post_file '/sleep_data/' subname '/session_' sess_format '/ses-' sess_format '_task-sleep_desc-preproc_Atlas_s2.55_regress.dtseries.nii'];
        if exist(cifti_file)
            sleep_sub_stage_cnt(stagei+1, subi) = 1 + sleep_sub_stage_cnt(stagei+1, subi);
            if isempty(sleep_sub_file_info{stagei+1, subi})
                sleep_sub_file_info{stagei+1, subi} = {cifti_file};
                sleep_sub_stage_info{stagei+1, subi} = {sess_format};
            else
                tai = sleep_sub_file_info{stagei+1, subi}; tai{end+1} = cifti_file;
                sleep_sub_file_info{stagei+1, subi} = tai;

                tai = sleep_sub_stage_info{stagei+1, subi}; tai{end+1} = sess_format;
                sleep_sub_stage_info{stagei+1, subi} = tai;
            end
        end

    end

end



EO_sub_ses_info = cell(1, length(sublist));
EO_sub_file_info = cell(1, length(sublist));
EO_sub_ses_cnt = zeros(1, length(sublist));

for subi = 1:length(sublist)
    subname = num2str(sublist(subi));
    for eoi = 1:2
        cifti_file = sprintf([EO_post_file '/sleep_data/%d/session_01/ses-01_task-sleep_desc-preproc_Atlas_s2.55_regress.dtseries.nii'], eoi, sublist(subi));
        
        if exist(cifti_file)
            EO_sub_ses_cnt(subi) = EO_sub_ses_cnt(subi) + 1; 
           if isempty(EO_sub_file_info{1, subi})
                EO_sub_file_info{1, subi} = {cifti_file};
                EO_sub_ses_info{1, subi} = {num2str(eoi)};
            else
                ta = EO_sub_file_info{1, subi};
                ta{end+1} = cifti_file;
                EO_sub_file_info{1, subi} = ta;

                ta = EO_sub_ses_info{1, subi};
                ta{end+1} = num2str(eoi);
                EO_sub_ses_info{1, subi} = ta;
            end 
        end

    end

end


sleep_sub_stage_info(6, :) = EO_sub_ses_info;
sleep_sub_file_info(6, :) = EO_sub_file_info;
sleep_sub_stage_cnt(6, :) = EO_sub_ses_cnt;

sleep_fd_info = cell(size(sleep_sub_stage_info));
sleep_dvars_info = cell(size(sleep_sub_stage_info));
% 
% sleep_fd_list = cell(6, 1);
% sleep_dvars_list = cell(6, 1);

for stai = 1:6
    for subi = 1:length(sublist)
        subid = sublist(subi);
        temp_fd = {};
        temp_dvars = {};
        for si = 1:length(sleep_sub_stage_info{stai, subi})
            sess_txt = sleep_sub_stage_info{stai, subi};
            if stai == 6
                regress_f = sprintf(EO_post_regress_file, sess_txt{si}, subid, sess_txt{si}, subid, sess_txt{si});
            else
                regress_f = sprintf(Sleep_300s_regress_file, subid, sess_txt{si}, subid, sess_txt{si});
            end
            regress_table = readtable(regress_f, 'FileType', 'text');
            fd = regress_table.framewise_displacement;
            dvars = regress_table.dvars;
            fd(isnan(fd)) = []; fd(fd == 0) = [];
            dvars(isnan(dvars)) = []; dvars(dvars == 0) = [];
            temp_fd{end+1} = fd;
            temp_dvars{end+1} = dvars;

%             sleep_fd_list{stai}(end+1) = fd;
%             sleep_dvars_list{stai}(end+1) = dvars;
        end
        sleep_fd_info{stai, subi} = temp_fd;
        sleep_dvars_info{stai, subi} = temp_dvars;
    end
end





save([ work_code_dir '/sleep_subinfo.mat'], 'sleep_sub_stage_info', 'sleep_sub_file_info', 'sleep_sub_stage_cnt', 'sublist', 'sleep_fd_info', 'sleep_dvars_info');





