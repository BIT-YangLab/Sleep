
function smooth_merge_sleep_cifti_post_dtseries(runs_sub_list, proj_dir, file_stem_t, start_s, end_s)

% runs_sub_list='/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/Gordon2016_mod/example/sublist_eo_all_new_test_post.txt';
% 
% proj_dir='/nd_disk3/guoyuan/sleep/sleep_eo_ciftify_2/postpreprocess_ciftify';




sublist = importdata(runs_sub_list);

% sleep_data/3089/session_01/ses-01_task-sleep_desc-preproc_Atlas_s0_regress.dtseries.nii
temp_surf_file = ['/nd_disk3/guoyuan/Jinlong/cifti_smooth_roi_' file_stem_t];

for subi = start_s:end_s%length(sublist)
    tic
    subname = num2str(sublist(subi));sprintf('%02d', subi);
    disp(['process ' num2str(subi) '/' num2str(length(sublist))]);
    fmri_list=[proj_dir '/scripts/rfMRI/subject_lists/allsub_MSM_ICA_FIX_fMRI_list/' subname '_rfMRI_list.txt'];
    if ~exist(fmri_list)
        continue
    end

    [sess, run, fmri, nii] = textread(fmri_list, '%s%s%s%s');
    data_dir = split(fmri{1}, 'Results');
    surf_dir = [data_dir{1} '/fsaverage_LR32k/'];

    for s_i = 1:length(sess)
        ses_n = sess{s_i};
        post_ret_cifti = [proj_dir '/sleep_data/' subname '/session_' ses_n '/' ses_n '_Atlas_s0_regress.dtseries.nii'];
        post_smooth_cifti = [proj_dir '/sleep_data/' subname '/session_' ses_n '/' ses_n '_Atlas_smooth_regress.dtseries.nii'];
        if ~exist(post_ret_cifti)
            continue
        end
        a1 = ft_read_cifti(post_ret_cifti);
        a1.dtseries(~isnan(a1.dtseries)) = 1;
%         ft_write_cifti(temp_surf_file, a1, 'parameter', 'dtseries');
        cmd = [ 'wb_command -cifti-smoothing ' post_ret_cifti ' 2.55 3.4 ' char("'") 'COLUMN' char("'") ' ' post_smooth_cifti ...
            '  -left-surface ' surf_dir 'sub-' subname '.L.midthickness.32k_fs_LR.surf.gii ' ...
            ' -right-surface ' surf_dir 'sub-' subname '.R.midthickness.32k_fs_LR.surf.gii ' ... 
            ' -cifti-roi ' temp_surf_file '.dtseries.nii ' ];
        disp(cmd);
        system(cmd);
    end
    toc
end

end