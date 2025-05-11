function new_label_assign = subcortex_parcellations_cleanup(orig_label, varargin)
% reference: Lisa et al., 2019
    p = inputParser;
    addParameter(p, 'cerebellum_voxel_idx', []);
    addParameter(p, 'dist_subcortex_subcortex', []);
    addParameter(p, 'neighbor_subcortex', []);
    parse(p, varargin{:});

    cerebellum_voxel_idx = p.Results.cerebellum_voxel_idx;
    dist_subcortex_subcortex = p.Results.dist_subcortex_subcortex;
    neighbor_subcortex = p.Results.neighbor_subcortex;

    template_dlabel = '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/resource/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_netassignments_LR.dlabel.nii';
    setenv('CBIG_CODE_DIR', '/nd_disk3/guoyuan/Jinlong/CBIG');
    addpath([ getenv('CBIG_CODE_DIR') '/utilities/matlab/fslr_matlab' ]);
    addpath([ getenv('CBIG_CODE_DIR') '/external_packages/SD/SDv1.5.1-svn593/BasicTools/' ]);
    lh_targ_mesh = CBIG_read_fslr_surface('lh', 'fs_LR_32k','inflated','medialwall.annot');
    rh_targ_mesh = CBIG_read_fslr_surface('rh', 'fs_LR_32k','inflated','medialwall.annot');
    dlabel_struct = ft_read_cifti(template_dlabel);

    position_whole = dlabel_struct.pos; 
    position_whole(1:32492, :) = lh_targ_mesh.vertices';
    position_whole(32493:64984, :) = rh_targ_mesh.vertices';

    if isempty(cerebellum_voxel_idx)
        MNI_ras2vox_file(1:96854, './test.nii.gz');
        nifti_base_info = MRIread('./test.nii.gz');
        system('rm ./test.nii.gz');
        template_MNI152_cortex = '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/resource/MNI152_cortex_mask.nii.gz';
        MNI152_struct = MRIread(template_MNI152_cortex);
        nzk = MNI152_struct.vol ~= 0 & nifti_base_info.vol ~= 0;
        overlap_idx = nifti_base_info.vol(nzk);
        cerebellum_voxel_idx = zeros(96854, 1);
        cerebellum_voxel_idx(ans) = overlap_idx;
    end 
    if isempty(dist_subcortex_subcortex)
        dist_subcortex_subcortex = compute_vertices_dist(position_whole(64985:end, :), position_whole(64985:end, :));
    end 
    if isempty(neighbor_subcortex)
        neighbor_subcortex = cell(96854-64984, 1);
        for voxel_i = 1:(96854-64984)
            neighbor_list = find_neighbor_voxel(voxel_i, dist_subcortex_subcortex, position_whole(64985:end, :));
            neighbor_subcortex{voxel_i} = neighbor_list;
        end
    end

    new_label_clean = orig_label;
    orig_exist_ali_idx = find(orig_label ~= 0);

    % step1: remove cerebellar voxels within 2mm of the cortex from the initial network assignment
    cerebellum_voxel_idx_flg = cerebellum_voxel_idx == 1;
    % cerebellar Left:10, Right:11  (fs_LR_32k)
    cerebellar_idx = (dlabel_struct.brainstructure == 10) + (dlabel_struct.brainstructure == 11);
    cerebellum_voxel_idx_flg(cerebellar_idx == 0) = 0; 
    new_label_clean(cerebellum_voxel_idx_flg ~= 0) = 0;

    % step2: remove isolate voxel and parcels:
    % removed isolated single voxel parcels that did not share a network assignment with any 
    % adjacent voxels, and parcels of size 2–4 voxels that did not have a counterpart
    % with the same network assignment within a 2mm radius in the contralateral hemisphere. (Ji et al., 2019.)
    % Parcels which shared a corner and had a continuous contralateral counterpart were combined
    
    neigbor_net_cluster_idx = {};
    cluster_count_idx = zeros(96854-64984, 1);
    nk_cluster = 0;
    for voxel_i = 1:(96854-64984)
        if orig_label(voxel_i+64984) == 0
            continue;
        end
        
        curr_net = orig_label(voxel_i+64984);
        neibor_net_list = orig_label(neighbor_subcortex{voxel_i}+64984);
        neibor_net_list(neibor_net_list < 1) = 0;

        % single voxel 
        if nnz(neibor_net_list == curr_net) == 0
            new_label_clean(voxel_i+64984) = 0;
            
            continue;
        end

        % multi-voxel 
        if sum(neibor_net_list == curr_net) > 0 
            if cluster_count_idx(voxel_i) ~= 0
                nk_cluster_m = cluster_count_idx(voxel_i);
            else
                
                nk_cluster_m = nk_cluster + 1;
            end
            neighbor_temp_cluster_idx = [voxel_i neighbor_subcortex{voxel_i}(neibor_net_list == curr_net)];
            neighbor_temp_cluster_idx = unique(neighbor_temp_cluster_idx);
            if nnz(cluster_count_idx(neighbor_temp_cluster_idx)) == 0
                neigbor_net_cluster_idx{nk_cluster_m} = neighbor_temp_cluster_idx;
                cluster_count_idx(neigbor_net_cluster_idx{nk_cluster_m}) = nk_cluster_m;
                nk_cluster = nk_cluster_m;
            else
                t1 = cluster_count_idx(neighbor_temp_cluster_idx); t1(t1 < 1) = [];
                cluster_count_idx(neighbor_temp_cluster_idx) = min(t1);
                neigbor_net_cluster_idx{min(t1)} = [neigbor_net_cluster_idx{min(t1)} neighbor_temp_cluster_idx];
                neigbor_net_cluster_idx{min(t1)} = unique(neigbor_net_cluster_idx{min(t1)});
            end
        end
    end

    % remove parcels with 2-4 voxels
    cluster_dis = zeros(100000, 1);
    for cluster_i = 1:length(neigbor_net_cluster_idx)
        parcels_voxel = neigbor_net_cluster_idx{cluster_i};
        cluster_dis(length(parcels_voxel)) = cluster_dis(length(parcels_voxel)) + 1;
        if length(parcels_voxel) < 5 
            new_label_clean(parcels_voxel+64984) = 0;
        end
    end

    % used nearest-neighbour interpolation to reassign the voxels removed from 
    % network assignment in the previous steps
    

    new_label_assign = new_label_clean;
    
    % reassign parcels with 2-4
    for cluster_i = 1:length(neigbor_net_cluster_idx)
        if length(neigbor_net_cluster_idx{cluster_i}) >= 5
            continue;
        end
        voxel_idx = neigbor_net_cluster_idx{cluster_i};
        neibor_net_list = [];
        for voxel_j = 1:length(voxel_idx)
            net1 = new_label_clean(neighbor_subcortex{voxel_idx(voxel_j)}+64984);
            net1(net1 < 1) = [];
            neibor_net_list = [neibor_net_list; net1];
        end
        if isempty(neibor_net_list)
            % some voxel has no neighbor node
            try_times = 10;
            while try_times
                voxel_2nd_neighbor = [];
                for voxel_j = 1:length(voxel_idx)
                    net1 = new_label_assign(neighbor_subcortex{voxel_idx(voxel_j)}+64984);
                    net1(net1 < 1) = [];
                    neibor_net_list = [neibor_net_list; net1];
                    voxel_2nd_neighbor = [ voxel_2nd_neighbor neighbor_subcortex{voxel_idx(voxel_j)}];
                    voxel_2nd_neighbor = unique(voxel_2nd_neighbor);
                end
                if isempty(neibor_net_list)
                    voxel_idx = voxel_2nd_neighbor;
                    try_times = try_times - 1;
                else
                    break;
                end
    
            end
            
        end
        if isempty(neibor_net_list)
            disp(['lost neighbor voxel: cluster ' num2str(cluster_i) ]);
            new_label_assign(voxel_idx+64984) = orig_label(voxel_idx+64984);
            continue;
        end
        new_label_assign(voxel_idx+64984) = mode(neibor_net_list);
    end

    % reassign single voxel and removed voxels in cerebellar
    new_interpolation_idx = intersect(orig_exist_ali_idx, find(new_label_assign == 0));
    for voxel_i = 1:length(new_interpolation_idx)
        voxel_idx = new_interpolation_idx(voxel_i)-64984;
        neibor_net_list = [];
        
        net1 = new_label_clean(neighbor_subcortex{voxel_idx}+64984);
        net1(net1 < 1) = [];
        neibor_net_list = [neibor_net_list; net1];
        
        if isempty(neibor_net_list)
            % some voxel has no neighbor node
            try_times = 10;
            while try_times
                voxel_2nd_neighbor = [];
                for voxel_j = 1:length(voxel_idx)
                    net1 = new_label_assign(neighbor_subcortex{voxel_idx(voxel_j)}+64984);
                    net1(net1 < 1) = [];
                    neibor_net_list = [neibor_net_list; net1];
                    voxel_2nd_neighbor = [ voxel_2nd_neighbor neighbor_subcortex{voxel_idx(voxel_j)}];
                    voxel_2nd_neighbor = unique(voxel_2nd_neighbor);
                end
                if isempty(neibor_net_list)
                    voxel_idx = voxel_2nd_neighbor;
                    try_times = try_times - 1;
                else
                    break;
                end
    
            end
            
        end
        if isempty(neibor_net_list)
            disp(['lost neighbor voxel: cluster ' num2str(cluster_i) ]);
            new_label_assign(voxel_idx+64984) = orig_label(voxel_idx+64984);
            continue;
        end
        new_label_assign(voxel_idx+64984) = mode(neibor_net_list);
    end
end

