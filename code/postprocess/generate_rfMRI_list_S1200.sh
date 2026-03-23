#!/bin/bash
proj_dir=$1;
proj_dir1=$1/out;          # hard coded; the folder storing your downloaded data
data_dir=$3
# origin ABCC subject_rsfMRI filepath

runs_sub_list=$2;    # hard coded

source /nd_disk3/guoyuan/sleep/a_bash/Gordon_method/postprocess/postprocess_config_test7t.txt

mkdir -p $proj_dir
mkdir -p $proj_dir/scripts/rfMRI/subject_lists/allsub_MSM_ICA_FIX_fMRI_list
cd $proj_dir

echo ${runs_sub_list}
# cat ${runs_sub_list}
# exit 0
#for subname in ??????; do
# session=1;
cat ${runs_sub_list} | while read subname; do
	echo ${subname}
	num_run=1
	rfmri_list=$proj_dir/scripts/rfMRI/subject_lists/allsub_MSM_ICA_FIX_fMRI_list/${subname}_rfMRI_list.txt
	echo $rfmri_list
	if [ -f $rfmri_list ]; then
		rm $rfmri_list
	fi
	
    echo ${data_dir}/sub-${subname}
	if [ ! -d ${data_dir}/sub-${subname} ]; then
		continue
	fi

	#sess_all_num=`ls ${data_dir}/sub-${subname}/MNINonLinear/Results | grep ${fmri_list_dir_grep} | wc -l`
    #sess_all_dir=(`ls ${data_dir}/sub-${subname}/MNINonLinear/Results`)

	echo "ls ${data_dir}/sub-${subname}/ | grep ${fmri_list_dir_grep} | wc -l"
	
	sess_all_num=`ls ${data_dir}/sub-${subname}/ | grep ${fmri_list_dir_grep} | wc -l`

    sess_all_dir=(`ls ${data_dir}/sub-${subname}/ | grep ${fmri_list_dir_grep}`)
    
	for session in $(seq 0 ${sess_all_num})
	do
		sess_number=`echo ${session} | awk '{printf("%02d\n",$0)}'`
        
        dir_1=${sess_all_dir[$session]}
        
		
        #file_cifti=`ls ${data_dir}/sub-${subname}/MNINonLinear/Results/$dir_1 | grep ${fmri_list_cifti_grep}`
        #file_nifti=`ls ${data_dir}/sub-${subname}/MNINonLinear/Results/$dir_1 | grep ${fmri_list_nifti_grep}`

		file_cifti=`ls ${data_dir}/sub-${subname}/$dir_1/func | grep ${fmri_list_cifti_grep}`
        file_nifti=`ls ${data_dir}/sub-${subname}/$dir_1/func | grep ${fmri_list_nifti_grep}`
		
        echo $file_cifti
		
		

		#ICA_fixed_fMRI=${data_dir}/sub-${subname}/MNINonLinear/Results/ses-${sess_number}_task-sleep_desc-preproc/ses-${sess_number}_task-sleep_desc-preproc_Atlas_s0.dtseries.nii
		#native_fMRI=${data_dir}/sub-${subname}/MNINonLinear/Results/ses-${sess_number}_task-sleep_desc-preproc/ses-${sess_number}_task-sleep_desc-preproc.nii.gz
		eval "ICA_fixed_fMRI=$fmri_list_cifti_format"
        eval "native_fMRI=$fmri_list_nifti_format"
        
        echo ==============================
        echo cifti: $ICA_fixed_fMRI
		echo nifti: $native_fMRI
		echo ==============================
        
		if [ $process_nifti -eq 1 ]; then
			if [ -f $ICA_fixed_fMRI ]; then
				
				echo "${dir_1} $num_run $native_fMRI $native_fMRI" >> $rfmri_list
				num_run=$(($num_run+1))
			fi
		else
			if [ -f $ICA_fixed_fMRI ]; then
				
				echo "${dir_1} $num_run $ICA_fixed_fMRI $native_fMRI" >> $rfmri_list
				num_run=$(($num_run+1))
			fi
		fi
	done
	

done
