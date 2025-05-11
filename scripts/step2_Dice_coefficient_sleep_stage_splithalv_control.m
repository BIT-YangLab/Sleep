% Author: Jinlong Li and Guoyuan Yang, BIT.
% compute Dice's coefficient for control group
% 

clc;clear;


work_dir = '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/';
ret_dir = '/nd_disk3/guoyuan/Jinlong/sleep_ret/HFR_ret/';

% network_label_all sleep_sub_file_info sleep_sub_stage_cnt sublist sleep_sub_stage_info  filter_sub_idx 
load('/nd_disk3/guoyuan/Jinlong/sleep_ret/HFR_ret/IndiPar_net_multi_7/Network_ind_subcortex_cleanup5min.mat');

load('/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/resource/fs_LR_32k_no_medial_wall_index.mat');  %no_medial_wall_index




load([ ret_dir '/subject_info.mat']);

min_subn = min(sum(sleep_sub_stage_cnt(1:4, :) ~= 0, 2));
min_subn = floor(min_subn / 2)*2;

template_dlabel = '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/resource/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_netassignments_LR.dlabel.nii';

dlabel_struct = ft_read_cifti(template_dlabel);
a = ft_read_cifti(template_dlabel);
cortex_mask_lh = a.brainstructure == 1;
cortex_mask_rh = a.brainstructure == 2;
cerebellum_mask_lh = a.brainstructure == 10;
cerebellum_mask_rh = a.brainstructure == 11;
thalamus_mask_lh = a.brainstructure == 20;
thalamus_mask_rh = a.brainstructure == 21;
striatum_mask_lh = a.brainstructure == 3 | a.brainstructure == 8 | a.brainstructure == 18;
striatum_mask_rh = a.brainstructure == 4 | a.brainstructure == 9 | a.brainstructure == 19;

mask_list = {cortex_mask_lh, cortex_mask_rh, cerebellum_mask_lh, cerebellum_mask_rh, thalamus_mask_lh, thalamus_mask_rh};
mask_name_list = {'cortex_lh', 'cortex_rh', 'cerebellum_lh', 'cerebellum_rh', 'thalamus_lh', 'thalamus_rh'};

mask_list_w = {cortex_mask_lh | cortex_mask_rh, cerebellum_mask_lh | cerebellum_mask_rh, thalamus_mask_lh | thalamus_mask_rh, striatum_mask_lh | striatum_mask_rh, ...
    cortex_mask_lh | cortex_mask_rh | cerebellum_mask_lh | cerebellum_mask_rh | thalamus_mask_lh | thalamus_mask_rh | striatum_mask_lh | striatum_mask_rh};
mask_name_list_w = {'cortex', 'cerebellum', 'thalamus', 'striatum', 'wholebrain'};

mask_idx = 1;

network_label_all = network_label_subcortex_lisa;

network_n = 7;
vertex_n = 96854;
random_time = 1000; % split half 

load([work_dir '/subject_control_select_7.mat']);

for random_i = 1:random_time
    for stai = 1:4

        % find non repetitive sublist
        sub_exist = control_sub_list_select{stai};
        try_times = 384200;
        while try_times
            sub1_idx = randperm(length(sub_exist), length(sub_exist)/2);
            sub1_list = sub_exist(sub1_idx);
            sub2_list = setdiff(sub_exist, sub1_list);
            unique_sublist_flag = 1;
            for rand_his_i = 1:random_i-1
                if ~isempty(history_random_list{stai, rand_his_i})
                    continue;
                end
                if nnz(sub1_list == history_random_list{stai, rand_his_i}) == length(sub1_list) || ...
                    nnz(sub2_list == history_random_list{stai, rand_his_i}) == length(sub2_list)
                    unique_sublist_flag = 0;
                    break;
                end
            end
            if unique_sublist_flag
                history_random_list{stai, random_i} = sub1_idx;
                break;
            else
                try_times = try_times - 1;
                if mod(try_times, 10000) == 0
                    disp(['process stage:' num2str(stai) ' time:' num2str(try_times)]);
                end
            end

            
        end
        if try_times == 0
            
            break;
        end

        temp_random_sublist = {sub1_list, sub2_list};
        temp_label_id = zeros(vertex_n, 2);
        for split_i = 1:2
            temp_label_sta = zeros(vertex_n, network_n);
            nc = 0;
            for sub_idx = 1:length(temp_random_sublist{split_i})
                subi = temp_random_sublist{split_i}(sub_idx);
                
                if isempty(network_label_all{stai, subi})
                    continue;
                end
                label_tmp = network_label_all{stai, subi}{1};
                for ni = 1:network_n
                    idx = find(label_tmp == ni);
                    temp_label_sta(idx, ni) = temp_label_sta(idx, ni) + 1;
                end
                nc = nc + 1;
            end
            temp_label_sta = temp_label_sta ./ nc;
            [~, colIndex] = max(temp_label_sta, [], 2);
            colIndex(setdiff(1:64984, no_medial_wall_index)) = 0;
            temp_label_id(:, split_i) = colIndex;
        end
        network_ret{stai, random_i} = temp_label_id;
    end
end


Dice_list_mean_cell = cell(4, 1);

for mask_idx = 1:5

    DICE_list_mean = zeros(4, 4);
    DICE_list_label = zeros(4, 4, network_n);

    total_valid_n = 0;
    for random_i = 1:random_time
        bk_flg = 1;
        for stai = 1:4
            if isempty(network_ret{stai, random_i})
                bk_flg = 0;
                break;
            end
        end
        if bk_flg == 1
            total_valid_n = total_valid_n + 1;
            network_ret1(:, total_valid_n) = network_ret(:, random_i);
        end
    end

    DICE_list_mean = zeros(4, 4);
    rk = zeros(4, 4);

    for random_i = 1:total_valid_n
        
        
        for stai = 1:4
            for staj = 1:4
                if isempty(network_ret1{stai, random_i}) || isempty(network_ret1{staj, random_i})
                    continue;
                end
                if stai == staj
                    [dice_list, label_list] = compute_dice_coef(network_ret1{stai, random_i}(mask_list_w{mask_idx}, 1), network_ret1{stai, random_i}(mask_list_w{mask_idx}, 2));
                else
                    [dice_list, label_list] = compute_dice_coef(network_ret1{stai, random_i}(mask_list_w{mask_idx}, 1), network_ret1{staj, random_i}(mask_list_w{mask_idx}, 1));
                end
                DICE_list_mean(stai, staj) = DICE_list_mean(stai, staj) + nanmean(dice_list);
                
                rk(stai, staj) = rk(stai, staj) + 1;
            end
        end
    end

    DICE_list_mean = DICE_list_mean ./ rk;
    % control
    DICE_list_mean = DICE_list_mean(1:4, 1:4);

    Dice_list_mean_cell{mask_idx} = DICE_list_mean;
end

% control
caxis_list = [
    0.65, 0.85;
    0.25, 0.5;
    0.10, 0.50;
    0.10, 0.50;
    0.50, 0.8;
];



for mask_idx = 1:4
    DICE_list_mean = Dice_list_mean_cell{mask_idx};
    draw_fig(DICE_list_mean, caxis_list(mask_idx, 1), caxis_list(mask_idx, 2), ['/nd_disk3/guoyuan/Jinlong/sleep_ret/HFR_ret/pic/DICE_control_heatmap_mean_net_' num2str(7)  '_' mask_name_list_w{mask_idx} '.png'])
    draw_fig_clean(DICE_list_mean, caxis_list(mask_idx, 1), caxis_list(mask_idx, 2), ['/nd_disk3/guoyuan/Jinlong/sleep_ret/HFR_ret/pic/DICE_control_heatmap_mean_net_' num2str(7)  '_clean_' mask_name_list_w{mask_idx} '.png'])
end

function draw_fig_clean(data, min_v, max_v, out_name)
    % draw fig
    stageNames = {'Awake', 'N1', 'N2', 'N3', 'REM'};

    pMatrix = data;
    % Set significance level
    alpha1 = 0.001;
    alpha2 = 0.01;
    alpha3 = 0.05;

    % Create a color depth map
    fig = figure;
    fig.Position = [100, 100, 1000, 800];
    imagesc(pMatrix); % Visualize using ImageSC

    
    colorMap = [
        0.933, 0.980, 0.918;
        0.9, 0.9, 0.8;    
        1, 0.7, 0.7;      
        1, 0, 0;          
        0.827, 0.110, 0.110;
        0.75, 0, 0;        
        
    ];

    % Subdivide color mapping to achieve gradient
    nSteps = 256;  % Add details and generate more color transitions
    gradients = linspace(0, 1, nSteps);  % Create a gradient value from 0 to 1

    % Interpolate color mapping
    smoothColorMap = [];
    for i = 1:size(colorMap, 1) - 1
        % Obtain the RGB values of each adjacent color pair
        startColor = colorMap(i, :);
        endColor = colorMap(i + 1, :);
        
        % Interpolate these two colors to generate a gradient color
        gradColors = [linspace(startColor(1), endColor(1), nSteps/size(colorMap, 1))', ...
                    linspace(startColor(2), endColor(2), nSteps/size(colorMap, 1))', ...
                    linspace(startColor(3), endColor(3), nSteps/size(colorMap, 1))'];
        
        % Add gradient colors to smooth color mapping
        smoothColorMap = [smoothColorMap; gradColors];
    end

    % Reset color mapping
    colormap(smoothColorMap);  

    caxis([min_v max_v]);     


    set(gca, 'FontSize', 18);       % 设置坐标轴刻度字体大小


    hold on;

    % Manually draw black borders for each grid
    for i = 1:size(pMatrix, 1)
        for j = 1:size(pMatrix, 2)
            rectangle('Position', [j-0.5, i-0.5, 1, 1], ...
                    'EdgeColor', 'k', 'LineWidth', 1); % Draw a black border
        end
    end
    hold off;

    

    set(gca, 'xtick', [], 'xticklabel', [])
    set(gca, 'ytick', [], 'yticklabel', [])
    export_fig(out_name, '-m6', '-q100');


end


function draw_fig(data, min_v, max_v, out_name)
    % draw fig
    stageNames = {'Awake', 'N1', 'N2', 'N3', 'REM'};

    pMatrix = data;
    % Set significance level
    alpha1 = 0.001;
    alpha2 = 0.01;
    alpha3 = 0.05;

    % Create a color depth map
    fig = figure;
    fig.Position = [100, 100, 1000, 800];
    imagesc(pMatrix); % Visualize using ImageSC

    
    colorMap = [
        0.933, 0.980, 0.918;
        0.9, 0.9, 0.8;    
        1, 0.7, 0.7;      
        1, 0, 0;          
        0.827, 0.110, 0.110;
        0.75, 0, 0;        
        
    ];

    % Subdivide color mapping to achieve gradient
    nSteps = 256;  % Add details and generate more color transitions
    gradients = linspace(0, 1, nSteps);  % Create a gradient value from 0 to 1

    % Interpolate color mapping
    smoothColorMap = [];
    for i = 1:size(colorMap, 1) - 1
        % Obtain the RGB values of each adjacent color pair
        startColor = colorMap(i, :);
        endColor = colorMap(i + 1, :);
        
        % Interpolate these two colors to generate a gradient color
        gradColors = [linspace(startColor(1), endColor(1), nSteps/size(colorMap, 1))', ...
                    linspace(startColor(2), endColor(2), nSteps/size(colorMap, 1))', ...
                    linspace(startColor(3), endColor(3), nSteps/size(colorMap, 1))'];
        
        % Add gradient colors to smooth color mapping
        smoothColorMap = [smoothColorMap; gradColors];
    end


    colormap(smoothColorMap);  % Set gradient color mapping

    caxis([min_v max_v]);     % Set color range
    title('Dice Coefficients');
    xlabel('sleep stage', 'FontSize', 18);
    ylabel('sleep stage', 'FontSize', 18);

    % Set the scale and labels for the x-axis and y-axis
    xticks(1:4);                     % Set the x-axis scale to 1 to 5
    yticks(1:4);                     % Set the y-axis scale to 1 to 5
    xticklabels(stageNames);         % Set x-axis scale label
    yticklabels(stageNames);         % Set y-axis scale label

    set(gca, 'FontSize', 18);       % Set the font size of the coordinate axis scale

    % Add grid and saliency annotations
    hold on;
    for i = 1:size(pMatrix, 1)
        for j = 1:size(pMatrix, 2)
            % Display numerical values in black font
            text(j, i, sprintf('%.3f', pMatrix(i, j)), ...
                'HorizontalAlignment', 'center', ...
                'Color', 'black', 'FontSize', 12, 'FontWeight', 'bold');
        end
    end
    % Manually draw black borders for each grid
    for i = 1:size(pMatrix, 1)
        for j = 1:size(pMatrix, 2)
            rectangle('Position', [j-0.5, i-0.5, 1, 1], ...
                    'EdgeColor', 'k', 'LineWidth', 1); % Draw a black border
        end
    end
    hold off;

    % Add color dial
    c = colorbar;  


    c.Ticks = [];  
    c.TickLabels = {};  


    caxis([min_v max_v]);  


    c.Ticks = [min_v, max_v]; 
    c.TickLabels = {num2str(min_v), num2str(max_v)};  


    c.Position = [0.91, 0.2, 0.03, 0.6];  

    set(gca, 'xtick', [], 'xticklabel', [])
    set(gca, 'ytick', [], 'yticklabel', [])
    export_fig(out_name, '-m6', '-q100');



end