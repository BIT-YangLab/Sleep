% Author: Jinlong Li and Guoyuan Yang, BIT.
% construct individual atlas for cerebral cortex by Wang2015(Wang et al., 2015)
% for each subject, construct atlas with a single session 

clc;clear;



work_dir = '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method';

flg_have_profile = 1;
flg_have_inter_variability = 1;
network_n = 7;
flg_have_grp_template = 1;



%%% add path of the subfolders
global ProgramPath
ProgramPath = [ work_dir '/HFR_ai_mod' ] ;
RESOURCE_DIR = [ProgramPath '/resource'];

addpath(genpath([ProgramPath '/code']));
addpath(genpath([ProgramPath '/scripts']));



dr_fs32k_fold = ''; % data reading file (fs_LR_32k)
dr_out = '/nd_disk3/guoyuan/Jinlong/sleep_ret/HFR_ret' ; % output dir
sub_list = ''; %subject list
confidence_threshold = 1.2;
numIter = 100;
combineLeftRight = 0; % combine hemispheres in iterative functional process

%% create information list
format_eo_cifti = '/nd_disk3/guoyuan/sleep/sleep_eo_ciftify_%d/postpreprocess_ciftify/sleep_data/%d/session_01/ses-01_task-sleep_desc-preproc_Atlas_s2.55_regress.dtseries.nii';
format_sleep_cifti = '/nd_disk3/guoyuan/sleep/sleep_300s_ciftify_post_new/postpreprocess_ciftify/sleep_data/%d/session_%s/ses-%s_task-sleep_desc-preproc_Atlas_s2.55_regress.dtseries.nii';

% sleep_sub_stage_info sleep_sub_file_info sleep_sub_stage_cnt sublist
load([ProgramPath '/sleep_subinfo.mat']);



%% create stable group template networks 
CBIG_CODE_DIR = '/nd_disk3/guoyuan/Jinlong/CBIG';
addpath([ CBIG_CODE_DIR, '/utilities/matlab/speedup_gradients']);
addpath([ CBIG_CODE_DIR, '/utilities/matlab/fslr_matlab']);
addpath([ CBIG_CODE_DIR, '/external_packages']);
addpath([ CBIG_CODE_DIR, '/external_packages/SD/SDv1.5.1-svn593/BasicTools/']);
addpath([ CBIG_CODE_DIR, '/utilities/matlab/speedup_gradients/utilities/']);
addpath([ CBIG_CODE_DIR, '/external_packages/SD/SDv1.5.1-svn593/kd_tree/']);
addpath([ CBIG_CODE_DIR, '/external_packages/matlab/default_packages/others/']);
addpath([ CBIG_CODE_DIR, '/utilities/matlab/FC/']);
addpath([ CBIG_CODE_DIR, '/utilities/matlab/stats/']);
addpath([ CBIG_CODE_DIR, '/utilities/matlab/surf/']);
addpath([ CBIG_CODE_DIR, '/external_packages/matlab/default_packages/graph_cut/matlab/']);
addpath([ CBIG_CODE_DIR, '/external_packages/matlab/default_packages/matlab_bgl/']);
addpath([ CBIG_CODE_DIR, '/external_packages/matlab/default_packages/mtimesx_20110223/']);
addpath([ CBIG_CODE_DIR, '/utilities/matlab/parcellation/']);
addpath([ CBIG_CODE_DIR, '/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Code/lib']);
addpath([ CBIG_CODE_DIR, '/external_packages/matlab/default_packages/DSP/']);
addpath([ CBIG_CODE_DIR, '/utilities/matlab/DSP/']);
addpath([ CBIG_CODE_DIR, '/utilities']);
addpath([ CBIG_CODE_DIR, '/utilities/matlab/predictive_models/PFM']);
addpath([ CBIG_CODE_DIR, 'stable_projects', 'predict_phenotypes', 'Kong2023_GradPar']);
addpath([ CBIG_CODE_DIR, '/utilities/matlab/predictive_models/utilities']);
addpath([ CBIG_CODE_DIR, '/utilities/matlab/utilities']);
addpath([ CBIG_CODE_DIR, '/utilities/matlab/predictive_models/KernelRidgeRegression']);
addpath([ CBIG_CODE_DIR, '/utilities/matlab/predictive_models/LinearRidgeRegression']);
addpath([ CBIG_CODE_DIR, '/utilities/matlab/stats']);

out_dir = [dr_out '/Group_template'];
seed_mesh = 'fs_LR_900';
targ_mesh = 'fs_LR_32k';
if ~exist(out_dir)
    mkdir([out_dir '/all_profile']);
end

load([RESOURCE_DIR '/fs_LR_32k_no_medial_wall_index.mat']);  %no_medial_wall_index
if network_n == 7
    ref_group = load([RESOURCE_DIR  '/HNU_example_clusters07_scrub_cifti_post_smooth_standard.mat']);
else
    ref_group = load([RESOURCE_DIR  '/HCP_40sub_1000iter_17nets_cen_sm4.mat']);
end
label_ref = [ref_group.lh_labels; ref_group.rh_labels];


temp_txt_profile = [dr_out '/test_profile.txt'];

if flg_have_profile == 0
    raw_i = 1;

    for subi = 1:length(sublist)
        
        disp(['compute RSFC profile: ' num2str(raw_i) '/' num2str(length(sublist)*2)]);
        subname = sublist(subi);
        if isempty(sleep_sub_file_info{6, subi})
            continue;
        end
        file_l = sleep_sub_file_info{6, subi};
        for ki = 1:length(file_l)
            cifti_file = file_l{ki};
            fid = fopen(temp_txt_profile, 'w');
            fprintf(fid, '%s\n', cifti_file);
            fclose(fid);
            profile_file = [ out_dir '/all_profile/' num2str(subname) '_sid' num2str(ki) '.mat' ];
            raw_i = raw_i + 1;
            if exist(profile_file)
                continue;
            end
    
            CBIG_ComputeCorrelationProfile(seed_mesh,targ_mesh, profile_file, 'NONE', '0.1', temp_txt_profile, 'NONE', 'NONE', '0');
        end
    end
    
end


% iterative select subjects to get group atlas

if flg_have_grp_template == 0
    n_iter_time = 10;
    group_atlas_label = zeros(64984, n_iter_time);
    new_sublist = repmat(sublist, 2, 1);

    wake_file_all = {};
    for subi = 1:length(sublist)
        subname = sublist(subi);
        if isempty(sleep_sub_file_info{6, subi})
            continue;
        end
        file_l = sleep_sub_file_info{6, subi};
        for ki = 1:length(file_l)
            profile_file = [ out_dir '/all_profile/' num2str(subname) '_sid' num2str(ki) '.mat' ];
            wake_file_all{end+1} = profile_file;
        end
    end

    rng(5904, "twister");
    for itn = 1:n_iter_time
        k_idx = randperm(length(wake_file_all), 150);
        temp_subfile_rand = wake_file_all(k_idx);
        num_data = 0;
        for kidx = 1:length(temp_subfile_rand)
            
            profile_file = temp_subfile_rand{kidx};
            if exist(profile_file)
                num_data = num_data + 1;
                
                profile_data = load(profile_file);

                if(num_data == 1)
                    avg_profile = profile_data;
                else
                    avg_profile.profile_mat = avg_profile.profile_mat + profile_data.profile_mat;
                end
            end 

        end
        profile_mat = avg_profile.profile_mat./num_data;
        save([out_dir '/avg_profile1.mat'],'profile_mat','-v7.3');
        CBIG_VonmisesSeriesClustering_fix_bessel_randnum_bsxfun(targ_mesh, 'NONE', network_n, [out_dir '/group_large1.mat'], ...
        [out_dir '/avg_profile1.mat'], 'NONE', 0, 1000, 0, 100, 1);

        

        a1 = load([out_dir '/group_large1.mat']);
        label = [a1.lh_labels; a1.rh_labels];

        [output, assign, cost, dice_overlap] = CBIG_HungarianClusterMatch(label_ref, label);
        group_atlas_label(:, itn) = output;

    end


    group_network_match = zeros(network_n, 64984);
    for itn = 1:n_iter_time
        for neti = 1:network_n
            idx = find(group_atlas_label(:, itn) == neti );
            group_network_match(neti, idx) = group_network_match(neti, idx) + 1;
        end
    end
    group_network_match = group_network_match ./ n_iter_time;
    save([ RESOURCE_DIR '/Template/Yeo' num2str(network_n) 'net_grp_template.mat'], "group_network_match")
else
    load([ RESOURCE_DIR '/Template/Yeo' num2str(network_n) 'net_grp_template.mat'])
end

if flg_have_inter_variability == 0
    %% create information list
    
    load([ProgramPath '/sleep_subinfo.mat']);

    % generate inter-subjects' variability connectome across all vertex
    profile_file_list = {};
    for subi = 1:length(sublist)
        subname = sublist(subi);
        if isempty(sleep_sub_file_info{6, subi})
            continue;
        end
        file_l = sleep_sub_file_info{6, subi};
        for ki = 1:length(file_l)
            profile_file = [ out_dir '/all_profile/' num2str(subname) '_sid' num2str(ki) '.mat' ];
            profile_file_list{end+1} = profile_file;
        end
    end

    variability = generate_inter_variability(profile_file_list, 10);
    save([ RESOURCE_DIR '/sleep_inter_individual_variability.mat'], "variability");
end



%% Wang2015 iterative functional region


vertex_n = 32492;


% Load distance matrix, which will be used for smooth the vertex-vertex
% correlation
fs32k_Firstadjacent_file = [RESOURCE_DIR '/fs32k_Firstadjacent_vertex.mat'];

% Load individual variability obtained from muller et al, Neuro, 2013
load([ RESOURCE_DIR '/sleep_inter_individual_variability.mat']); % variability
Prior_Variability = variability';%Prior_Variability:0.5~0.8
Prior_SNR = ones(1, vertex_n*2);
GRP_template = [ RESOURCE_DIR '/Template/Yeo' num2str(network_n) 'net_grp_template.mat' ];

network_label_all = cell(size(sleep_sub_stage_info));

% stagei == 6 actually means eo1 datasets
for stagei = 1:6
    OutPath = [dr_out '/IndiPar_net_multi' num2str(network_n) '/st_' num2str(stagei) ];
    if ~exist(OutPath)
        mkdir(OutPath);
    end
    for s = 1:length(sublist)
        sub = sublist(s);
        disp(['Individual Parcellation ' num2str(s) ':' num2str(sub) '\n']);
        outdir = [OutPath '/' num2str(sub)];
        if exist([outdir '/Network_ind_ret.mat'])
            load([outdir '/Network_ind_ret.mat']);
            network_label_all{stagei, s} = network_cur_iter;
            continue;
        end
        temp_sleep_info = sleep_sub_file_info{stagei, s};
        temp_sleep_stage = sleep_sub_stage_info{stagei, s};

        if max(size(temp_sleep_info)) == 0
            continue;
        end

        network_cur_iter_list = {};
        for sesi = 1:length(temp_sleep_info)
            cifti_f = temp_sleep_info{sesi};
            sesn = temp_sleep_stage{sesi};
            if ~exist(cifti_f)
                continue;
            end
            cifti_data = ft_read_cifti(cifti_f);
            raw_BOLD_set = cifti_data.dtseries;
            if isempty(raw_BOLD_set)
                continue; 
            end
            network_cur_iter = Func_iterativeParcellation_mod(raw_BOLD_set, Prior_Variability, Prior_SNR, outdir, numIter, ...
            GRP_template, network_n, confidence_threshold, fs32k_Firstadjacent_file, combineLeftRight);
            network_cur_iter_list{end+1} = network_cur_iter;
        end
        network_label_all{stagei, s} = network_cur_iter_list;
    end
end

save([OutPath '/../Network_ind_ret.mat'], "network_label_all", "sleep_sub_file_info", "sublist", "sleep_sub_stage_cnt", "sleep_sub_stage_info");
















