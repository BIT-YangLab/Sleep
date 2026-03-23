
# /nd_disk3/guoyuan/Jinlong/postproc/submit_HCP_S1200_postprocessing_with_censoring_without_interpolation.sh 
# .. /nd_disk3/guoyuan/sleep/a_bash/Gordon_method/Gordon2016_mod/example/sublist_eo_all_new_test_post.txt 
# ../nd_disk3/guoyuan/sleep/sleep_eo_ciftify_2/postpreprocess_ciftify_nogsr 
# ../nd_disk3/guoyuan/sleep/sleep_eo_ciftify_2/out/ciftify/
# runs_sub_list="/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/Gordon2016_mod/example/sublist_eo_all_new_test_post.txt";
# proj_dir="/nd_disk3/guoyuan/sleep/sleep_eo_ciftify_2/postpreprocess_ciftify";
# data_dir="/nd_disk3/guoyuan/sleep/sleep_eo_ciftify_2/out/ciftify/"

source /nd_disk3/guoyuan/sleep/a_bash/Gordon_method/postprocess/postprocess_config_test7t.txt

runs_sub_list=${sub_list_file}
proj_dir=${postprocess_out_dir}
data_dir=${preprocess_data_ciftify_dir}

sdir_path="${proj_dir}/sleep_data"
FD_ths=${FD_threshold}
DVARS_ths=${DVARS_threshold}
sh_file="${postprocess_work_dir}/HCP_S1200_postprocessing_general_with_censor_without_interpolation.csh"

cd ${postprocess_work_dir}

${postprocess_work_dir}/generate_rfMRI_list_S1200.sh $proj_dir $runs_sub_list $data_dir

# no need to run script, existing dvars' file 
#/nd_disk3/guoyuan/Jinlong/postproc/submit_HCP_S1200_generate_DVARS.sh $runs_sub_list $proj_dir  $data_dir

p_cnt=0

for subname in `cat $runs_sub_list`; do
	date
	echo $subname
	
	# if [ -f ${sdir_path}/${subname}/session_1/ses-01_task-sleep_desc-preproc_Atlas_s0_regress.dtseries.nii ]; then
	# 	continue
	# else
	# 	echo "do $subname"
	# fi


	log=${proj_dir}/scripts/logs_S1200/$subname.log;
	if [ -f "$log" ] ; then
		rm -rf ${log}
	fi
	

	mkdir -p ${proj_dir}/scripts/logs_S1200

	
	fmri_list="${proj_dir}/scripts/rfMRI/subject_lists/allsub_MSM_ICA_FIX_fMRI_list/${subname}_rfMRI_list.txt"
    echo $fmri_list
	if [ ! -f $fmri_list ]; then
		echo "skip sub-$subname"
		continue
	fi

	cmd="csh ${sh_file} -s $subname -s-dir ${sdir_path} -input-type MSM_ICA_FIX -fmri-list ${fmri_list} -project_name MSM_reg_wbsgrayordinatecortex -WBS graycortex -FD-threshold ${FD_ths} -DV-threshold ${DVARS_ths} -discard-segment-length 5 -discard-run 75 >& ${log}"
	# eval $cmd
	# exit 0
	if [ $p_cnt -eq 3 ]; then
		eval $cmd 
		p_cnt=0
	else
		cmd_n="${cmd} &"
		eval $cmd_n 
		p_cnt=$((p_cnt+1))
	fi
	date

done




