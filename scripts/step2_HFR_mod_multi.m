% Author: Jinlong Li and Guoyuan Yang, BIT.
% construct individual atlas for cerebral cortex by Wang2015(Wang et al., 2015)
% for each subject, construct atlas with all sessions

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
load([dr_out '/subject_info.mat']);
load([dr_out '/sublist_filter_new.mat']);



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

%% Yeo2011 group profile
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


%% iterative select subjects to get group atlas of Yeo 7network
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


%% create inter variability Mueller el al., 2013
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
load([RESOURCE_DIR '/fs32k_Firstadjacent_vertex.mat']);

% Load individual variability obtained from muller et al, Neuro, 2013
load([ RESOURCE_DIR '/sleep_inter_individual_variability.mat']); % variability
Prior_Variability = variability';%Prior_Variability:0.5~0.8
Prior_SNR = ones(1, vertex_n*2);

network_label_all = cell(size(sleep_sub_stage_info));

% stagei == 6 actually means eo1 datasets
for stagei = 1:5
    OutPath = [dr_out '/IndiPar_net' num2str(network_n) '_plus/st_' num2str(stagei) ];
    if ~exist(OutPath)
        mkdir(OutPath);
    end
    for s = 1:length(sublist)
        sub = sublist(s);
        disp(['Individual Parcellation ' num2str(s) ':' num2str(sub) '\n']);
        if exist([OutPath '/' num2str(sub) '/Network_ind_ret.mat'])
            load([OutPath '/' num2str(sub) '/Network_ind_ret.mat']);
            network_label_all{stagei, s} = network_cur_iter;
            continue;
        end
        if ~sleep_sub_filter_cnt(stagei, s)
            continue;
        end
        temp_sleep_info = sleep_sub_file_info{stagei, s};

        if max(size(temp_sleep_info)) == 0
            continue;
        end

        raw_BOLD_set = [];
        for sesi = 1:length(temp_sleep_info)
            cifti_f = temp_sleep_info{sesi};
            if ~exist(cifti_f)
                continue;
            end
            cifti_data = ft_read_cifti(cifti_f);
            % contral
            %raw_BOLD = cifti_data.dtseries(:, 11:160);
            raw_BOLD = cifti_data.dtseries;
            vol=bsxfun(@minus,raw_BOLD,mean(raw_BOLD,2)); % row wise demeaning
            raw_BOLD=bsxfun(@rdivide,vol,std(vol',1)'); % row wise standardization 
            raw_BOLD_set = [raw_BOLD_set raw_BOLD];
        end
        if isempty(raw_BOLD_set)
            continue; 
        end
        lhData = raw_BOLD_set(1:32492, :);
        rhData = raw_BOLD_set(32493:64984, :);

        seedDatalh = [];
        seedDatarh = [];

        n = size(lhData,1);
        m = size(rhData,1);

        %%%%%%%%%% Think about the priors, Could include variability and SNR
        if range(Prior_Variability) > 0
            Prior_Variability = 0.4 + 0.6*(Prior_Variability- min(Prior_Variability))/(max(Prior_Variability) - min(Prior_Variability)); % normalize the range to 0.4 ~1. Therefore the inv will be between 1~2.5.
        end

        var_lh = Prior_Variability(1:n);
        var_rh = Prior_Variability(n+1:end);

        varInv_lh = 1./var_lh;
        varInv_rh = 1./var_rh;

        if range(Prior_SNR) >0
            Prior_SNR = 0.4+  0.6*(Prior_SNR- min(Prior_SNR))/(max(Prior_SNR) - min(Prior_SNR)) ; % normalize the range to 0.4 ~1. Therefore the inv will be between 1~2.5.
        end

        SNR_lh = Prior_SNR(1:n);
        SNR_rh = Prior_SNR(n+1:end);

        % ---------------------------------------------------
        %% Iterative parcellation
        % ---------------------------------------------------
        
        mkdir([OutPath '/' num2str(sub)]);
        for cnt = 1:numIter
            disp(['compute individual atlas: ' num2str(sub) '  iter:' num2str(cnt)]);tic;
            mkdir([OutPath '/' num2str(sub) '/Iter_' num2str(cnt)]);


            if cnt==1

                % group_network_match
                load([ RESOURCE_DIR '/Template/Yeo' num2str(network_n) 'net_grp_template.mat' ]);
                group_network_match_lh = group_network_match(:, 1:32492);
                group_network_match_rh = group_network_match(:, 32493:end);

                vol = group_network_match_lh(1, :);
                ventLh = find(vol>0);
                GrpNetlh{1}= ventLh;
                vol = group_network_match_rh(1, :);
                ventRh = find(vol>0);
                GrpNetrh{1}= ventRh;
                for i2=1:network_n  % get the seed waveforms based on Thomas' parcellation, and weight it by inv(Variability)
                    vol = group_network_match_lh(i2+1, :);
                    idx =find(vol>0);

                
                    seedDatalh(i2,:)= varInv_lh(idx)*lhData(idx,:); % weight the group map using the inverse of individual difference
                    GrpNetlh{i2+1} = idx;

                    vol = group_network_match_rh(i2+1, :);
                    idx =find(vol>0);

                
                    seedDatarh(i2,:)= varInv_rh(idx)*rhData(idx,:); % weight the group map using the inverse of individual difference
                    GrpNetrh{i2+1} = idx;
        
                end
                GrpSeedDatalh = seedDatalh;
                GrpSeedDatarh = seedDatarh;

            else
                load([OutPath '/' num2str(sub) '/Iter_' num2str(cnt-1) '/Iterative_group_network_match.mat' ]);


                for i2 = 1:network_n  % get the seed waveforms based on the last parcellation
                    vol = group_network_match_lh(i2+1, :);
                    vol(isnan(vol)) = 0;
                    idx = find(vol>=confidence_threshold);
                    if isempty(idx)
                        maxx = max(max(max(vol)));
                        idx = find(vol==maxx);
                    end
                    seedDatalh(i2,:) = SNR_lh(idx)*lhData(idx,:); % weight the individual signal based on SNR

                    vol = group_network_match_rh(i2+1, :);
                    vol(isnan(vol)) = 0;
                    idx = find(vol>=confidence_threshold);
                    if isempty(idx)
                        maxx = max(max(max(vol)));
                        idx = find(vol==maxx);
                    end

                    seedDatarh(i2,:)= SNR_rh(idx)*rhData(idx,:);

                end
            end

            %% Weight in the group seed in each iteration, should throw in individual variability map as weight in the future

            if cnt>1
                seedDatalh = seedDatalh + GrpSeedDatalh/(cnt-1);
                seedDatarh = seedDatarh + GrpSeedDatarh/(cnt-1);
            end


            % Combine the same network of left hemi and right hemi?,
            % If not, uncomment the following 3 lines

            % if combine the two hemisphere, the impact of other hemisphere is decreasing during iteration
            if (combineLeftRight)
                tmp = seedDatalh;
                seedDatalh = seedDatalh+seedDatarh/(cnt+2);
                seedDatarh = seedDatarh+tmp/(cnt+2);
    
            end



            % compute vertex to seed correlation for all vertices

            cValuelh = zeros(size(lhData,1),size(seedDatalh,1));
            cValuerh = zeros(size(rhData,1),size(seedDatarh,1));


            data = [seedDatalh;lhData];
            tmp = corrcoef(data');
            cValuelh = tmp(1:size(seedDatalh,1),(end-vertex_n+1):end)'; % 2562*seeds

            data = [seedDatarh;rhData];
            tmp = corrcoef(data');
            cValuerh = tmp(1:size(seedDatarh,1),(end-vertex_n+1):end)'; % 2562*seeds 

            cValuelh = 0.5*log((1+cValuelh)./(1-cValuelh));
            cValuerh = 0.5*log((1+cValuerh)./(1-cValuerh));

            cValuelh(isnan(cValuelh(:)))=0;
            cValuerh(isnan(cValuerh(:)))=0;



            % Smooth, decrease the noise
            for i = 1:vertex_n
                cValuelhS(i,:) = mean(cValuelh(fs32k_Firstadjacent_vertex_lh{i},:),1);
                cValuerhS(i,:) = mean(cValuerh(fs32k_Firstadjacent_vertex_rh{i},:),1);

            end

            cValuelh = cValuelhS;
            cValuerh = cValuerhS;

        
            % Further weight in the group map * inv(Variability) by adding correlation coefficient of 0~ 0.5 according to inv(Variability).

            for i = 1:network_n
                idx = GrpNetlh{i+1};
                cValuelh(idx, i) = cValuelh(idx, i) + (((varInv_lh(idx)-1)/3)/cnt)';

                idx = GrpNetrh{i+1};
                cValuerh(idx, i) = cValuerh(idx, i) + (((varInv_rh(idx)-1)/3)/cnt)';
            end

            % --------------------------------------------------------
            %%  Determine the network membership of each vertex
            % ------------  Left hemisphere---------------------------
        
            data=cValuelh(:,1:network_n);
            parc_membership = zeros(size(data, 1), 1);
            parc_confidence = zeros(size(data, 1), 1);
            for v=1:size(data,1)

                [cor idx] = sort(data(v,:),'descend');
                parc_membership(v) = idx(1);
                parc_confidence(v) = cor(1)/cor(2);

            end

            network_s = zeros(size(data, 1), network_n);
            network_cur_iter = zeros(vertex_n*2, 1);
        
            for n =1:network_n
                network = 0*parc_membership;
                confid = 0*network;
                network(find(parc_membership==n))= 1;
                confid(find(parc_membership==n))= parc_confidence(find(parc_membership==n));

                network(ventLh) = 0; % mask out the ventrical and useless areas in the midline
                network(isnan(network)) = 0;
                network_s(:, n) = network;

                group_network_match_lh(n+1, :) = network.*confid;
            end
            save([OutPath '/' num2str(sub) '/Iter_' num2str(cnt) '/Network_ind_lh.mat'], "network_s");
            network_cur_iter(1:vertex_n) = compute_network_label(network_s, network_n);
            
            
            % --------------------------------------------------------
            %%  Determine the network membership of each vertex
            % ------------  Right hemisphere---------------------------

            data=cValuerh(:,1:network_n);
            parc_membership = zeros(size(data, 1), 1);
            parc_confidence = zeros(size(data, 1), 1);
            for v=1:size(data,1)

                [cor idx] = sort(data(v,:),'descend');
                parc_membership(v) = idx(1);
                parc_confidence(v) = cor(1)/cor(2);


            end

            network_s = zeros(size(data, 1), network_n);

            for n =1:network_n
                network = 0*parc_membership;
                confid = 0*network;
                network(find(parc_membership==n))= 1;
                confid(find(parc_membership==n))= parc_confidence(find(parc_membership==n));
                network(ventRh) = 0; % mask out the ventrical and useless areas in the midline

                network_s(:, n) = network;

                group_network_match_rh(n+1, :) = network.*confid;
            end

            save([OutPath '/' num2str(sub) '/Iter_' num2str(cnt) '/Network_ind_rh.mat'], "network_s");
            save([OutPath '/' num2str(sub) '/Iter_' num2str(cnt) '/Iterative_group_network_match.mat' ], "group_network_match_lh", "group_network_match_rh");
            network_cur_iter(1+vertex_n:end) = compute_network_label(network_s, network_n);

            if cnt > 1
                network_pre_iter = zeros(vertex_n*2, 1);
                load([OutPath '/' num2str(sub) '/Iter_' num2str(cnt-1) '/Network_ind_lh.mat']);
                network_pre_iter(1:vertex_n) = compute_network_label(network_s, network_n);
                load([OutPath '/' num2str(sub) '/Iter_' num2str(cnt-1) '/Network_ind_rh.mat']);
                network_pre_iter(1+vertex_n:end) = compute_network_label(network_s, network_n);
                percent_overlap = nnz(network_cur_iter == network_pre_iter) / length(network_cur_iter);
                if percent_overlap >= 0.98 || cnt == numIter
                    network_label_all{stagei, s} = network_cur_iter;
                    save([OutPath '/' num2str(sub) '/Network_ind_ret.mat'], "network_cur_iter");
                    break;
                end
            end

        toc
        end
    end
end

save([OutPath '/../Network_ind_ret.mat'], "network_label_all", "sleep_info", "sublist", "sub_stage_cell");

function network_all = compute_network_label(network_s, network_n)
    network_all = zeros(size(network_s, 1), 1);
    for ki = 1:network_n
        idx = find(network_s(:, ki) == 1);
        network_all(idx) = ki; 
    end
end














