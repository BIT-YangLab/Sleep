#!/bin/bash


## description
# This is the sample code for the Sleep project, 
# outlining the entire processing and invocation 
# workflow along with descriptive comments. 
# It is divided into three parts: 
# first, Data Processing, where raw data undergoes pre- and post-processing, 
# and individual sleep stages are identified and 
# normalized based on EEG staging info. 
# Second, the Project Workflow, which calls the interfaces 
# for each module (refer to specific scripts and docs for logic). 
# Finally, Data and Visuals, featuring figures generated via 
# GraphPad, Connectome Workbench and plot_fig_subcortex(https://github.com/wd-veloce/plot_fig_subcortex),
# manually typeset in PowerPoint. Brief data descriptions are provided here; 
# please refer to the documentation for further details.

WORK_DIR=../

###########################################################
##                   data preprocess                     ##
###########################################################


# install docker-images: tigrlab/fmriprep_ciftify latest

# out_dir
mkdir -p /nd_disk3/guoyuan/sleep/sleep_eo_ciftify_2
mkdir -p /nd_disk3/guoyuan/sleep/sleep_eo_ciftify_2/tmp12

# running
docker run -ti --rm \
    -v /md_disk4/guoyuan/Sleep/Sleep_raw/Sleep_eo_2/:/data:ro \
    -v /nd_disk3/guoyuan/sleep/sleep_eo_ciftify_2:/out \
    -v $HOME/license.txt:/opt/freesurfer/license.txt \
    -v /nd_disk3/guoyuan/sleep/sleep_eo_ciftify_2/tmp12:/tmp \
    -w /tmp \
    tigrlab/fmriprep_ciftify:latest \
    /data /out/out \
    participant \
    --participant_label=3001,3002 \
    --n_cpus 8 \
    --ignore-fieldmaps \
    --fs-license /opt/freesurfer/license.txt \
    --fmriprep-args="--use-aroma --ignore-aroma-denoising-errors --fs-license /opt/freesurfer/license.txt"

###########################################################
##                   data post-process                   ##
###########################################################
# more details please see: ${WORK_DIR}/code/postprocess/
${WORK_DIR}/code/postprocess/submit_HCP_S1200_postprocessing_with_censoring_without_interpolation.sh

###########################################################
##       Participant data organization                   ##
###########################################################
# Given that the imaging data had already been divided into sessions 
# according to EEG-defined sleep stages before further processing, 
# the session data of each individual were subsequently organized 
# and summarized based on these EEG session annotations, 
# and the relevant statistics were calculated.

${WORK_DIR}/scripts/step1_check_subjects.m
${WORK_DIR}/scripts/step1_prepare_sleep_subject_info_analysis.m
${WORK_DIR}/scripts/step1_select_control_sublist.m

###########################################################
## Construct individualized whole-brain functional atlas ##
###########################################################

# Multi-session data from each individual were concatenated into a single run, and a brain atlas was obtained.
${WORK_DIR}/scripts/step2_HFR_mod_multi.m
${WORK_DIR}/scripts/step2_cleanup_subcortex_multi.m

# The code example for the control group is provided; the remaining procedures are similar.
${WORK_DIR}/scripts/step2_subcotex_mapping_multi_control.m

# A brain atlas was derived from each individual session of each participant.
${WORK_DIR}/scripts/step2_HFR_mod_single_session.m
${WORK_DIR}/scripts/step2_cleanup_subcortex.m

###########################################################
##      Validate the stability of the brain atlas        ##
###########################################################
${WORK_DIR}/scripts/step2_Dice_coefficient_sleep_stage_splithalv.m
${WORK_DIR}/scripts/step2_Dice_coefficient_sleep_stage_splithalv_control.m

# Figure S2, NMI
${WORK_DIR}/scripts/step3_compute_atlas_nmi.m

###########################################################
##          Group-average atlas visualization            ##
###########################################################

# For whole-brain atlas visualization, the results were generated in 
# both dlabel and NIfTI formats, and the network area of the group-average atlas was calculated.
# Fig. 1, Figure S1, Figure S2
${WORK_DIR}/scripts/step3_visualization_subcortex_atlas.m
${WORK_DIR}/scripts/step3_visualization_subcortex_atlas_control.m

###########################################################
##          network size and significance                ##
###########################################################

# Figure S3
${WORK_DIR}/scripts/step4_compute_network_size_lme_whole.m
# Figure S4
${WORK_DIR}/scripts/step4_compute_network_size_lme_whole_control.m

###########################################################
##             Sleep stage classification                ##
###########################################################
# Multi-class SVM classification using network topography as features. Fig. 3
${WORK_DIR}/scripts/step5_classify_sleep_stage_onehot_multiclass_7_whole.m

# Multi-class SVM classification using network size as features. Figure S10
${WORK_DIR}/scripts/step5_classify_sleep_stage_multiclass_7_whole.m

# Binary-class SVM classification using network topography as features. Figure S9
${WORK_DIR}/scripts/step5_classify_sleep_stage_onehot_binary_7_whole.m

###########################################################
##             network encroaching                       ##
###########################################################

# Fig. 2, Figure S5
${WORK_DIR}/scripts/step6_compute_network_encroaching.m

# We will use the code from Neurosynth. For details, 
# please refer to https://github.com/neurosynth/neurosynth
# example: decoder = decode.Decoder(dataset, method='roi')
# result = decoder.decode(['MNI_neurosynth_encroach.nii.gz'], save='decoding_results.txt')
${WORK_DIR}/scripts/step6_create_encroach_neurosynth_mni.m


###########################################################
##                    Graph analysis                     ##
###########################################################
# Figure S6, spring-embeded visualization in Gephi (https://gephi.org/)
${WORK_DIR}/scripts/step6_compute_network_spring_embeded.m
${WORK_DIR}/scripts/step6_compute_network_modularity.m



###########################################################
##             sleep quality and SWA                     ##
###########################################################
# Fig. 4
${WORK_DIR}/scripts/step7_swa_related_with_surface_area_wholebrain.m
${WORK_DIR}/scripts/step8_behavior_related_with_surface_area_whole_brain.m






