
References:
If you are using this repo please cite:
xxx

Usage:
The code utilized in this study include eight main steps:
1. prepare data and check information
step1_check_subjects.m: 
    Check subjects' information include subject id, session, file path
step1_prepare_sleep_subject_info_analysis.m: 
    Organize relevant data of individuals for analysis: demographic information, FD, swa, behavior measures
step1_select_control_sublist.m: 
    Select subjects as control group for validation
step1_compute_cortex_snr_weight.m:
    Compute snr map of cortex for mapping subcortex atlas

2. construct individual and group average atlas 
step2_HFR_mod_multi.m/step2_HFR_mod_single_session.m: 
    Construct individual atlas in cortex using Wang et al., 2015.
step2_subcortex_mapping.m/step2_subcortex_mapping_multi_all.m/step2_subcortex_mapping_multi_control.m:
    Mapping individual atlas in subcortex using Lisa et al., 2019.
step2_cleanup_subcortex.m/step2_cleanup_subcortex_multi.m:
    Clean up individual atlas in sucortex
step2_Dice_coefficient_sleep_stage_splithalv_control.m/step2_Dice_coefficient_sleep_stage_splithalv.m:
    Compute Dice's coefficient for group average atlas

3. visualization of group average atlas
step3_visualization_subcortex_atlas.m/step3_visualization_subcortex_atlas_control.m:
    Visualization of group_average atlas

4. compute network size
step4_compute_network_size_lme_whole.m/step4_compute_network_size_lme_whole_control.m:
    Calculate the network size of individual atlas and save the data in LME model input file format
compute_linear_mixed_model_significance.R:
    LME model 

5. predict sleep stage using network topography
step5_classify_sleep_stage_binary_whole.m:
    SVM binary classification
step5_classify_sleep_stage_multiclass_7_whole.m/step5_classify_sleep_stage_multiclass_7_control_whole.m:
    SVM multi classification
step5_classify_sleep_stage_multiclass_7_whole_hemi.m:
    SVM multi classification, using network topology features with the left and right hemispheres 
step5_classify_sleep_stage_RF_multiclass_7_whole.m:
    RF multi classification
step5_classify_sleep_stage_KNN_multiclass_7_whole.m:
    KNN multi classification
step5_classify_sleep_stage_ANN_multiclass_7_whole.m:
    ANN multi classification

6. compute encroaching 
step6_compute_network_encroaching.m:
    Calculate the encroaching regions and types of AMN network
step6_compute_network_encroaching_rsfc.m:
    Calculate the RSFC between the encroaching regions of the AMN network in each sleep stage compared to the awake stage
step6_compute_network_spring_embeded.m:
    Calculate the data for spring embeded plotting
step6_compute_network_modularity.m:
    Calculate the modularity

7. calculate the correlation between network area and SWA
step7_swa_related_with_surface_area_wholebrain.m:
    Calculate the correlation between the network area and SWA for all sleep stages compared to the awake stage
step7_swa_related_with_surface_area_single.m:
    Calculate the correlation between the network area and SWA during a single sleep stages compared to the awake stages

8. calculate the correlation between network area and PSQI 
step8_behavior_related_with_surface_area_whole_brain.m
    Calculate the correlation between the network area and PSQI during a single sleep stages compared to the awake stages



