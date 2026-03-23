function CBIG_HCP_generate_regressors_clean(subname, out_dir, project_name, sess, runt)
% Written by Guoyuan Yang

fid = fopen('/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/postprocess/postprocess_config_test7t.txt');
C = textscan(fid, '%s %s', 'Delimiter', '=', 'CommentStyle', '#');
fclose(fid);

keys = C{1};
vals = C{2};

for i = 1:length(keys)
    cfg.(keys{i}) = vals{i};
end

sess_f = strsplit(sess, cfg.sess_split); sess_f = sess_f{1};
confounds_all=sprintf(cfg.fmriprep_confounds_format, subname, sess_f, subname, sess_f);


confounds_data = readtable(confounds_all, 'FileType', 'text');

transx = confounds_data.trans_x;
transy = confounds_data.trans_y;
transz = confounds_data.trans_z;
rotx = confounds_data.rot_x;
roty = confounds_data.rot_y;
rotz = confounds_data.rot_z;
csf = confounds_data.csf;
wm = confounds_data.white_matter;
wbs = confounds_data.global_signal;


csf_d = [0; diff(csf)];
wm_d = [0; diff(wm)];
wbs_d = [0; diff(wbs)];
transx_d = [0; diff(transx)];
transy_d = [0; diff(transy)];
transz_d = [0; diff(transz)];
rotx_d = [0; diff(rotx)];
roty_d = [0; diff(roty)];
rotz_d = [0; diff(rotz)];


regressor_concat=[csf wm wbs transx transy transz rotx roty rotz csf_d wm_d wbs_d transx_d transy_d transz_d rotx_d roty_d rotz_d];

dlmwrite(fullfile(out_dir, subname,['session_', sess], 'scripts', ['rfMRI_', sess, '_run-', runt, '_regressor_concat.txt']), regressor_concat);
