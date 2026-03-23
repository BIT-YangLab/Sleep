function CBIG_HCP_S1200_motion_outlier_detection_sleep(subname,data_dir,project_name,sess,runt,FD_th,DV_th,discard_segment_length)
% Written by Ru(by) Kong


FD_th = str2double(FD_th);
DV_th = str2double(DV_th);
if( ~isempty(discard_segment_length) )
    discard = 1;
    discard_seg = str2double(discard_segment_length);
end

fid = fopen('/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/postprocess/postprocess_config_test7t.txt');
C = textscan(fid, '%s %s', 'Delimiter', '=', 'CommentStyle', '#');
fclose(fid);

keys = C{1};
vals = C{2};

for i = 1:length(keys)
    cfg.(keys{i}) = vals{i};
end

sess_f = strsplit(sess, cfg.sess_split); sess_f = sess_f{1}; disp(sess_f)
DV_FD_file=sprintf(cfg.fmriprep_confounds_format, subname, sess_f, subname, sess_f);

FD_DV=readtable(DV_FD_file, 'FileType', 'text');

curr_RelRMS = FD_DV.framewise_displacement(2:end);
% curr_RelRMS=str2num(char(curr_RelRMS));
curr_RelRMS = [0;curr_RelRMS];
DV_vect=FD_DV.dvars(2:end);
% DV_vect=str2num(char(DV_vect));
DV_vect = [0;DV_vect];


FD_vect = (curr_RelRMS <= FD_th);
DV_vect=(DV_vect<=DV_th);
DV_censor_vect = (DV_vect & [DV_vect(2:end,1); 1] & [1; DV_vect(1:end-1)] & [1;1; DV_vect(1:end-2)]);
FD_censor_vect = (FD_vect & [FD_vect(2:end,1); 1] & [1; FD_vect(1:end-1)] & [1;1; FD_vect(1:end-2)]);

DV_FD_censor_vect = [DV_censor_vect & FD_censor_vect];
if (discard == 1)
    DV_FD_censor_segment = diff([0; DV_FD_censor_vect; 0]);
    seg_start = find(DV_FD_censor_segment == 1);
    seg_end = find(DV_FD_censor_segment == -1)-1;

    for idx = 1:length(seg_start)
        if(seg_end(idx)-seg_start(idx) + 1 < discard_seg)
            DV_FD_censor_vect(seg_start(idx):seg_end(idx)) = 0;
        end
    end
end



dlmwrite(fullfile(data_dir, subname, ['session_', sess], 'scripts', ['rfMRI_', sess, '_run-', runt, '_FD', num2str(FD_th), '_DV', num2str(DV_th), '_censoring.txt']), DV_FD_censor_vect); % 

end