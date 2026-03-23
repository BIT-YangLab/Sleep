#! /bin/csh -f
##############################
## Written by Ru(by) Kong
##############################

set flag_WM = 0
set flag_CSF = 0
set flag_WBS = 0
set flag_MT = 0
set fMRI_list = ""
set glm_output_list = ""
set regressor_list = "" 
set discard_segment_length = 5
set DV_th = 75
set FD_th = 0.2
goto parse_args;
parse_args_return:

goto check_params;
check_params_return:
###############################
# check if matlab exists
###############################

# modify
#set MATLAB=`which $CBIG_MATLAB_DIR/bin/matlab`
# if ($status) then
# 	echo "ERROR: could not find matlab"
# 	exit 1;
# endif


foreach line (`grep -v '^#' /nd_disk3/guoyuan/sleep/a_bash/Gordon_method/postprocess/postprocess_config_test7t.txt`)
    set key = `echo $line | cut -d= -f1`

	if ($key == special_for_bash) then
		
		break
	endif
    set val = `echo $line | cut -d= -f2-`
	
    set $key = $val
end



set output_path=${subject_dir}/../;

# paths in if-else are hard coded; if here is changed, line 49-52 need to be changed accordingly
if ($input_type == MSM_ICA_FIX) then
	set fMRI_list = ${output_path}/scripts/rfMRI/subject_lists/MSM_ICA_FIX_regression_lists/fMRI_lists/${subject}_fmri_list.txt
	set glm_output_list = ${output_path}/scripts/rfMRI/subject_lists/MSM_ICA_FIX_regression_lists/output_lists/${subject}_output_list.txt
	set regressor_list = ${output_path}/scripts/rfMRI/subject_lists/MSM_ICA_FIX_regression_lists/regressor_lists/${subject}_reg_list.txt
	set censor_list = ${output_path}/scripts/rfMRI/subject_lists/MSM_ICA_FIX_regression_lists/censor_lists/${subject}_cen_list.txt

else if ($input_type ==  Sulc_ICA_FIX) then
	set fMRI_list = /md_disk4/guoyuan/LSRW_data/REST_post/scripts/rfMRI/subject_lists/regression_lists/fMRI_lists/${subject}_fmri_list.txt
	set glm_output_list = /md_disk3/guoyuan/CHCP_rest/scripts/rfMRI/subject_lists/regression_lists/output_lists/${subject}_output_list.txt
	set regressor_list = /md_disk3/guoyuan/CHCP_rest/scripts/rfMRI/subject_lists/regression_lists/regressor_lists/${subject}_reg_list.txt
	set censor_list = /md_disk3/guoyuan/CHCP_rest/scripts/rfMRI/subject_lists/regression_lists/regressor_lists/${subject}_cen_list.txt
else 
	echo "unable to recognize! input_type can only be MSM_ICA_FIX or Sulc_ICA_FIX"
	exit 1
endif

#set censor_list = ""

# following 4 lines correspond to line 31-44; hard coded
mkdir -p ${output_path}/scripts/rfMRI/subject_lists/MSM_ICA_FIX_regression_lists/fMRI_lists
mkdir -p ${output_path}/scripts/rfMRI/subject_lists/MSM_ICA_FIX_regression_lists/output_lists
mkdir -p ${output_path}/scripts/rfMRI/subject_lists/MSM_ICA_FIX_regression_lists/regressor_lists
mkdir -p ${output_path}/scripts/rfMRI/subject_lists/MSM_ICA_FIX_regression_lists/censor_lists

if(-e $fMRI_list) then
	rm -rf $fMRI_list
endif
if(-e $glm_output_list) then
	rm -rf $glm_output_list
endif
if(-e $regressor_list) then
	rm -rf $regressor_list
endif
if(-e $censor_list) then
	rm -rf $censor_list
endif

mkdir -p ${output_path}/scripts/rfMRI/subject_lists/$project_folder    # hard coded
set fmri_file = ${output_path}/scripts/rfMRI/subject_lists/$project_folder/${subject}_rfMRI_list.txt    # hard coded; 
rm -rf $fmri_file
cp -rf $fmri_orig_file $fmri_file

foreach fmri_run ("`cat $fmri_file`") 
	echo " Current fMRI file info: $fmri_run"
	set curr_sess = (`echo $fmri_run | awk '{printf $1}'`)
	set curr_runt = (`echo $fmri_run | awk '{printf $2}'`)
	set curr_fmri = (`echo $fmri_run | awk '{printf $3}'`)
	
	echo "[INFO]: Current session is $curr_sess"
	echo "[INFO]: Current runtimes is $curr_runt"
	echo "[INFO]: Current fMRI file is $curr_fmri"
	
	
	
	set project_output = $subject_dir/$subject/session_${curr_sess}
	mkdir -p $project_output
	set scripts_path = $project_output/scripts
	if(-e $scripts_path && $curr_runt == '1' ) then
		rm -rf $scripts_path
	endif
	mkdir -p $scripts_path
	
	#################################
	### generate motion outliers
	#################################
	echo "####################################################"
	echo "## start generate motion outliers!"
	echo "####################################################"
	#CBIG_HCP_motion_outlier_detection(data_dir,project_name,subname,sess,run,FD_th,DV_th,discard_segment_length)

	set output_censor_file = $scripts_path/rfMRI_${curr_sess}_run-${curr_runt}_FD${FD_th}_DV${DV_th}_censoring.txt

	if(-e $output_censor_file) then
			echo "${output_censor_file} regressors already generated!"
	else
		

		set cmd = ( matlab -nodesktop -nodisplay -nosplash -r '"' 'addpath(genpath('"'${postprocess_work_dir}'"'))'; CBIG_HCP_S1200_motion_outlier_detection_sleep $subject $subject_dir $project_folder $curr_sess $curr_runt $FD_th $DV_th $discard_segment_length; exit; '"' ); # hard coded
		
		
		echo "$cmd"
		eval $cmd
	endif
	if(! -e $output_censor_file) then
		echo "$output_censor_file"
		echo "Fail!"
		exit 1;
	endif
	echo $output_censor_file >> $censor_list
	cat $censor_list
	echo "####################################################"
	echo "## discard run with more than $rm_run_th% censored frames!"
	echo "####################################################"

	if ($rm_run == 1) then 
		echo "========= check if the number of outliers exceeds $rm_run_th% of total number of frames for each run ========="
		set outlier_file = $output_censor_file
		echo "[motion outlier detection]: outlier_file = $outlier_file"
		set num_zeros = `fgrep -o 0 ${outlier_file} | wc -l`
		echo "[motion outlier detection]: num_zeros = $num_zeros" 
		set num_ones = `fgrep -o 1 ${outlier_file} | wc -l`
		echo "[motion outlier detection]: num_ones = $num_ones"
		@ num_frames = $num_zeros + $num_ones
		echo "[motion outlier detection]: num_frames = $num_frames"
		@ prop = $num_zeros * 100 / $num_frames
		echo "[motion outlier detection]: prop = $prop"
	
		if( `echo "$prop > $rm_run_th" | bc` ) then
			echo "[RUN DISCARD] session $curr_sess run has more than $rm_run_th% outliers, remove this run from fMRI list."
			sed -i "/${curr_sess}_/d" $fmri_file
			sed -i "/_${curr_sess}_/d" $censor_list
		else
			echo "[RUN KEEP] session $curr_sess run has less than $rm_run_th% outliers, nothing change to this run."
		endif

		echo "====================== $rm_run_th% outliers threshold check finished ======================"
		echo ""
	endif	
	cat $censor_list
end



####################################################
## generate regressors and do glm regression
####################################################

# first check if the censor list is empty, if it is, then this subject will be skipped
set check_remain_run = (`cat $fmri_file`)
set check_remain_run_no_space = (`echo "${check_remain_run}" | tr -d '[[:space:]]'`)
if ( "$check_remain_run_no_space" == "" ) then
	echo "[SUBJECT DISCARD] $subject has no run left after censoring." 
	echo "HCP postprocessing finished!"
	exit 1
else
	foreach fmri_run ("`cat $fmri_file`") 
		echo " Current fMRI file info: $fmri_run"
		set curr_sess = (`echo $fmri_run | awk '{printf $1}'`)
		set curr_runt = (`echo $fmri_run | awk '{printf $2}'`)
		set curr_fmri = (`echo $fmri_run | awk '{printf $3}'`)
	
		echo "[INFO]: Current session is $curr_sess"
		echo "[INFO]: Current runtimes is $curr_runt"
		echo "[INFO]: Current fMRI file is $curr_fmri"
	
		set project_output = $subject_dir/$subject/session_${curr_sess}
		set scripts_path = $project_output/scripts

	
		#################################
		### generate regressors
		#################################
		set temp_filename=$scripts_path/rfMRI_${curr_sess}_run-${curr_runt}_regressor_concat.txt
		if(-e $temp_filename) then
			echo "${temp_filename} regressors already generated!"
		else
			echo "####################################################"
			echo "## start generate regressors!"
			echo "####################################################"


			set cmd = ( matlab -nodesktop -nodisplay -nosplash -r '"' 'addpath(genpath('"'${postprocess_work_dir}'"'))'; CBIG_HCP_generate_regressors_clean $subject $subject_dir $project_folder $curr_sess $curr_runt; exit; '"' ); # hard coded echo "$cmd"
			
			echo $cmd
			eval $cmd
		endif
		if(! -e $temp_filename) then
			echo "Fail!"
			exit 1;
		endif


		echo $curr_fmri >> $fMRI_list
		echo "-------info---"
		cat $censor_list

		set output_file = (`basename $curr_fmri | awk -F'.' '{printf $1}'`)
		if ( "$process_nifti" == "1") then
			set output_path = $project_output/${output_file}_regress.nii.gz
		else
			set output_path = $project_output/${output_file}_regress
		endif

		echo $output_path >> $glm_output_list
	
		set regressor_filename = $temp_filename
		echo $regressor_filename >> $regressor_list
	
		echo "####################################################"


	end

	echo "####################################################"
	echo "## start glm regression with censoring!"
	echo "####################################################"
	set poly = 1
	set per_run = 1
    
    echo fmri:$fMRI_list
    echo glm_out:$glm_output_list
    echo censor:$censor_list

	echo "addpath('$CBIG_dir'); CBIG_glm_regress_vol '$fMRI_list' '$glm_output_list' '$regressor_list' '$poly' '$censor_list' '$per_run'; exit;"


	matlab -nodesktop -nodisplay -nosplash -r "addpath('$CBIG_dir'); CBIG_glm_regress_vol '$fMRI_list' '$glm_output_list' '$regressor_list' '$poly' '$censor_list' '$per_run'; exit;"


	##########################################
	# Check output
	##########################################
	set fail_flag = 0
	foreach fmri_run ("`cat $fmri_file`") 
		echo " Current fMRI file info: $fmri_run"
		set curr_sess = (`echo $fmri_run | awk '{printf $1}'`)
		set curr_runt = (`echo $fmri_run | awk '{printf $2}'`)
		set curr_fmri = (`echo $fmri_run | awk '{printf $3}'`)
	
		echo "[INFO]: Current session is $curr_sess"
		echo "[INFO]: Current runtimes is $curr_runt"
		echo "[INFO]: Current fMRI file is $curr_fmri"
		echo "[CHECK]: start check output"
		set project_output = $subject_dir/$subject/session_${curr_sess}
		set check_output = `ls $project_output/*${check_file_stem}`
		echo "$project_output/*${check_file_stem}"
		if($#check_output == 0) then
			set fail_flag = 1
		endif
	end


	if ($fail_flag == 0) then
		echo "HCP postprocessing finished!"
	else
		echo "Fail!"
	endif
	exit 1;
endif

##########################################
# Parse Arguments 
##########################################	

parse_args:
set cmdline = "$argv";
while( $#argv != 0 )
	set flag = $argv[1]; shift;
	
	switch($flag)
		#subject name
		case "-s":
			if ( $#argv == 0 ) goto arg1err;
			set subject = $argv[1]; shift;
			breaksw	
		
		case "-s-dir":
			if ( $#argv == 0 ) goto arg1err;
			set subject_dir = $argv[1]; shift;
			breaksw	
		
		case "-project_name":
			if ( $#argv == 0 ) goto arg1err;
			set project_folder = $argv[1]; shift;
			breaksw	
			
		case "-fmri-list":    # the list you generated from "generate_rfMRI_list_S1200.sh"
			if ( $#argv == 0 ) goto arg1err;
			set fmri_orig_file = $argv[1]; shift;
			breaksw
		
	    case "-input-type":
			if ( $#argv == 0 ) goto arg1err;
			set input_type = $argv[1]; shift;
			breaksw
			
		case "-FD-threshold":
			if ( $#argv == 0 ) goto arg1err;
			set FD_th = $argv[1]; shift;
			breaksw	
		
		case "-DV-threshold":
			if ( $#argv == 0 ) goto arg1err;
			set DV_th = $argv[1]; shift;
			breaksw
		
		case "-discard-segment-length":
			if ( $#argv == 0 ) goto arg1err;
			set discard_segment_length = $argv[1]; shift;
			breaksw
		
		case "-discard-run":
			set rm_run = 1;
			if ( $#argv == 0 ) goto arg1err;
			set rm_run_th = $argv[1]; shift;
			breaksw
						
		case "-WM":
			set flag_WM = 1;
			breaksw
				
		case "-CSF":
			set flag_CSF = 1; 
			breaksw
			
		case "-WBS"
			set flag_WBS = 1;
			if ( $#argv == 0 ) goto arg1err;
			set WBS_type = $argv[1]; shift;
			breaksw
			
		case "-MT"
			set flag_MT = 1;
			breaksw		
			
		default:
			echo ERROR: Flag $flag unrecognized.
			echo $cmdline
			exit 1
			breaksw
	endsw
end
goto parse_args_return;


##########################################
# Check Parameters
##########################################

check_params:
if ( "$subject" == "" ) then
	echo "ERROR: subject not specified"
	exit 1;
endif
if ( "$subject_dir" == "" ) then
	echo "ERROR: path to subject folder not specified"
	exit 1;
endif
goto check_params_return;


##########################################
# ERROR message
##########################################	

arg1err:
  echo "ERROR: flag $flag requires one argument"
  exit 1

arg2err:
  echo "ERROR: flag $flag requires two arguments"
  exit 1

	




