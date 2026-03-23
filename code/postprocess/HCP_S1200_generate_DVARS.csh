#!/bin/csh
## Written by Ru(by) Kong

set output_dir = $2/scripts/dvars_S1200   # hard coded; any intermediate folder
set data_dir = $3   # hard coded; the folder you store downloaded data
set postproc_dir = $2   # hard coded; the folder with postprocessed data

#set subject_list = "/share/users/imganalysis/yeolab/data/HCP/S1200/individuals/scripts/rfMRI/subject_lists/subject_runnum_list/all_new_subject_list_sort.txt"
#set subjects = `cat $subject_list`
#set subjects = "995174"

#foreach subname ($subjects)
	set subname = $1
	# set curr_sess = 1
	# set curr_run = 1

	set fmri_file = ${postproc_dir}/scripts/rfMRI/subject_lists/allsub_MSM_ICA_FIX_fMRI_list/${subname}_rfMRI_list.txt    # hard coded; corresponds to the input of HCP_S1200_postprocessing_general_with_censor_without_interpolation.csh -fmri-list
	foreach fmri_run ("`cat $fmri_file`") 
		echo " Current fMRI file info: $fmri_run"
		set curr_sess = (`echo $fmri_run | awk '{printf $1}'`)
		set curr_run = (`echo $fmri_run | awk '{printf $2}'`)
		set curr_fmri = (`echo $fmri_run | awk '{printf $3}'`)
		set curr_fmri_native = (`echo $fmri_run | awk '{printf $4}'`)
		
		echo "[INFO]: Current session is $curr_sess"
		echo "[INFO]: Current run is $curr_run"
	
		mkdir -p $output_dir/${subname}/rfMRI_${curr_sess}_${curr_run}/tmp
		
		if(-e $output_dir/${subname}/rfMRI_${curr_sess}_${curr_run}/${subname}_dvars) then
			echo "DVARS is already generated!"
		else
			fsl_motion_outliers -i ${curr_fmri_native} -o $output_dir/${subname}/rfMRI_${curr_sess}_${curr_run}/${subname}_dvars_confound -s $output_dir/${subname}/rfMRI_${curr_sess}_${curr_run}/${subname}_dvars -p $output_dir/${subname}/rfMRI_${curr_sess}_${curr_run}/${subname}_dvars_fig -t $output_dir/${subname}/rfMRI_${curr_sess}_${curr_run}/tmp --dvars --nomoco
		endif
		
		#break
		mkdir -p $postproc_dir/sleep_data/$subname/QC
		cp $output_dir/${subname}/rfMRI_${curr_sess}_${curr_run}/${subname}_dvars $postproc_dir/sleep_data/$subname/QC/DVARS.txt
	end 
#end
