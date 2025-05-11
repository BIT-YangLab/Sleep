% Author: Jinlong Li and Guoyuan Yang, BIT.
% Fit the correlation curve between the percentage of network size 
% and swa for a single sleep stage of a single network
% 
clc;clear;

work_dir = '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/';
ret_dir = '/nd_disk3/guoyuan/Jinlong/sleep_ret/HFR_ret/';

% network_label_all sleep_info sublist sub_stage_cell 
load([ret_dir '/IndiPar_net_multi_7/Network_ind_ret8min.mat']);
load([ work_dir '/resource/fs_LR_32k_no_medial_wall_index.mat']);  %no_medial_wall_index
load([ work_dir '/sleep_subinfo.mat']);
% line 6: eo1
sleep_sub_stage_cnt(:, ~filter_sub_idx) = [];
sleep_sub_file_info(:, ~filter_sub_idx) = [];
sleep_sub_stage_info(:, ~filter_sub_idx) = [];
sublist(~filter_sub_idx) = [];
network_label_all(:, ~filter_sub_idx) = [];


load([ ret_dir '/subject_info.mat']);
addpath(genpath('/nd_disk3/guoyuan/Xinyu/plot_fig_subcortex'));
addpath('/nd_disk3/guoyuan/Xinyu/software/cifti-matlab-master');
addpath('/nd_disk3/guoyuan/Xinyu/software/spm12');

load([ret_dir '/IndiPar_net_multi_7/Network_ind_subcortex_cleanup8min.mat']);


template_dlabel = '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/resource/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_netassignments_LR.dlabel.nii';

dlabel_struct = ft_read_cifti(template_dlabel);

cortex_mask_lh = dlabel_struct.brainstructure == 1;
cortex_mask_rh = dlabel_struct.brainstructure == 2;
cerebellum_mask_lh = dlabel_struct.brainstructure == 10;
cerebellum_mask_rh = dlabel_struct.brainstructure == 11;
thalamus_mask_lh = dlabel_struct.brainstructure == 20;
thalamus_mask_rh = dlabel_struct.brainstructure == 21;
striatum_mask_lh = dlabel_struct.brainstructure == 3 | dlabel_struct.brainstructure == 8 | dlabel_struct.brainstructure == 18;
striatum_mask_rh = dlabel_struct.brainstructure == 4 | dlabel_struct.brainstructure == 9 | dlabel_struct.brainstructure == 19;


mask_list = {cortex_mask_lh, cortex_mask_rh, cerebellum_mask_lh, cerebellum_mask_rh, thalamus_mask_lh, thalamus_mask_rh, striatum_mask_lh, striatum_mask_rh};
mask_name_list = {'cortex_lh', 'cortex_rh', 'cerebellum_lh', 'cerebellum_rh', 'thalamus_lh', 'thalamus_rh', 'striatum_lh', 'striatum_rh'};
mask_list_w = {cortex_mask_lh | cortex_mask_rh, cerebellum_mask_lh | cerebellum_mask_rh, thalamus_mask_lh | thalamus_mask_rh, striatum_mask_lh | striatum_mask_rh};
mask_name_list_w = {'cortex', 'cerebellum', 'thalamus', 'striatum'};

%% select different brain region
mask_idx = 2;


% line 6: eo1

network_n = 7;
vertex_n = 64984;

network_size_ret = cell(size(sleep_sub_stage_info));

network_label_all = network_label_subcortex_lisa;

for stai = 1:4
    
    nc = 0;
    for subi = 1:size(network_label_all, 2)
        if isempty(network_label_all{stai, subi})
            continue;
        end
        network_tmp_l = network_label_all{stai, subi};
        temp_size_sta = zeros(length(network_tmp_l), network_n);
        for sesi = 1:length(network_tmp_l)
            label_tmp = network_tmp_l{sesi};
            label_tmp = label_tmp(mask_list_w{mask_idx});
            for ni = 1:network_n
                idx = find(label_tmp == ni);
                if isempty(idx)
                    disp([num2str(stai) ' ' num2str(subi)])
                end
                temp_size_sta(sesi, ni) = length(idx);
            end
        end
        network_size_ret{stai, subi} = temp_size_sta./nnz(mask_list_w{mask_idx});
    end
end





network_size_ret1 = cell(4, size(network_size_ret, 2));
for subi = 1:size(network_size_ret, 2)
    if isempty(network_size_ret{1, subi})
        for st = 1:4
            network_size_ret1{st, subi} = [];
        end
        continue;
    end
    for st = 1:3
        if isempty(network_size_ret{st+1, subi})
            network_size_ret1{st, subi} = [];
            continue;
        end
        network_size_ret1{st, subi} = network_size_ret{st+1, subi} - network_size_ret{1, subi}(1, :);
    end
end


if network_n == 7
    label_name = {'Visual','Somatomotor','Dorsal attention', 'Action mode', 'Limbic', 'Frontoparietal', 'Default'};
elseif network_n == 17
    label_name = {'Auditory','Dorsal Attention A','Control A', 'Somatomotor A', 'Salience/VenAttn B', ...
'Default B', 'Default C', 'Visual C', 'Visual A', 'Dorsal Attention B', 'Temoral Pariental', 'Control B','Visual B',...
'Control C','Default A','Salience/VenAttn A', 'Somatomotor B'};
end

stageNames = {'Awake', 'N1', 'N2', 'N3', 'REM'};

network_size_variability = cell(network_n, 3);

R2_list = zeros(network_n, 4);



p_list_org = zeros(network_n, 3);
for networki = 1:network_n
    stage_size_variability = cell(4, 1);
    
    
    for st = 1:3
        kidx = 0;
        temp_size_variability = zeros(10, 3);
        for subi = 1:size(network_size_ret, 2)
            if isempty(network_size_ret1{st, subi})
                continue;
            end
            network_size_tmp = network_size_ret1{st, subi};
            for sesi = 1:size(network_size_tmp, 1)
                if isnan(swa_cell{st+1, subi}(sesi))
                    continue;
                end
                kidx = kidx + 1;
                temp_size_variability(kidx, :) = [network_size_tmp(sesi, networki) swa_cell{st+1, subi}(sesi) st+1];
            end
        end
        
        [r, p] = corr(temp_size_variability(:, 2), temp_size_variability(:, 1));
        p_list_org(networki, st) = p;
        network_size_variability{networki, st} = temp_size_variability;
    end
    
    
end


x_lim_list = {[0.15 0.6], [0.2 0.8], [0.25 0.85]};

for stagei = 1:3
    for networki = 4:4%network_n
        
        
        temp_size_variability = network_size_variability{networki, stagei};
        draw_plot(temp_size_variability, {'SWA'}, label_name{networki}, stageNames{1+1}, ['/nd_disk3/guoyuan/Jinlong/sleep_ret/HFR_ret/pic/swa_related_one_stage_' num2str(network_n) '_net' num2str(networki) '_' mask_name_list_w{mask_idx} '_' stageNames{stagei+1} '.png'], p_list_org(networki, stagei), x_lim_list{stagei});
        draw_plot_clean(temp_size_variability, {'SWA'}, label_name{networki}, stageNames{1+1}, ['/nd_disk3/guoyuan/Jinlong/sleep_ret/HFR_ret/pic/swa_related_one_stage_' num2str(network_n) '_net' num2str(networki) '_' mask_name_list_w{mask_idx} '_' stageNames{stagei+1} '_clean.png'], p_list_org(networki, 1), x_lim_list{stagei});
        
        
    end
end


    


function draw_plot(data, behavior_name, network_name, stage_name, out_name, p_fdr, x_lim)
    
    % Assuming the data is stored in a matrix data: the first column is Y, the second column is X, and the third column is the category
    Y = data(:, 1);  
    X = data(:, 2);  
    categories = data(:, 3);  

    % Calculate Pearson correlation coefficient r and P-value
    [r, p] = corr(X, Y);

    % Perform linear regression fitting
    model = fitlm(X, Y);

    % Extract the slope and intercept of the fit
    coeffs = model.Coefficients.Estimate;  

    % The fitted linear equation
    x_fit = linspace(min(X), max(X), 100);  
    [y_fit, ci] = predict(model, x_fit');  

    % Specify category color
    custom_colors = [0, 0, 0; 67, 205, 128; 238, 99, 99] / 255;  

    % Draw a scatter plot
    fig = figure;
    fig.Position = [100, 100, 800, 600];

    % Set the entire graphic background to white or transparent
    %fig.Color = 'white';  
    fig.Color = 'none'; 

    hold on;
    unique_categories = unique(categories);
    for i = 1:length(unique_categories)
        % Filter the points of the current category
        current_category = unique_categories(i);
        idx = categories == current_category;
        
        % Draw scatter points for the current category
        scatter(X(idx), Y(idx), 40, 'MarkerFaceColor', custom_colors(i, :), ...
            'MarkerEdgeColor', 'k');
    end

    % Draw confidence intervals
    fill([x_fit, fliplr(x_fit)], [ci(:, 1)', fliplr(ci(:, 2)')], ...
        [211, 211, 211] / 255, 'EdgeColor', 'none', 'FaceAlpha', 0.5);  

    % Draw a fitted line
    plot(x_fit, y_fit, 'r-', 'LineWidth', 2);  

    hold off;

    % Add labels
    xlabel(behavior_name, 'FontSize', 16);
    ylabel('variation of surface area', 'FontSize', 16);
    title(['Linear fit: ' network_name ], 'FontSize', 16);

    % Set coordinate axis
    set(gca, 'FontSize', 16);  
    set(gca, 'XColor', 'k', 'YColor', 'k', 'Box', 'off');  

    % Only retain the left Y-axis and bottom X-axis
    ax = gca;
    ax.LineWidth = 1.5;  
    ax.TickDir = 'in';   
    ax.XAxis.TickDirection = 'in'; 
    ax.YAxis.TickDirection = 'in';  

    % Mark the values of r and P in the bottom right corner
    text(max(X), min(Y), sprintf('r = %.2f\nP = %.3g', r, p_fdr), ...
        'FontSize', 16, 'Color', 'black', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
    xticks(0:0.2:1);
    yticks(-1:0.2:1);

    xlim(x_lim);


    export_fig(out_name, '-m6', '-q400');


end

function draw_plot_clean(data, behavior_name, network_name, stage_name, out_name, p_fdr, x_lim)
    
    % Assuming the data is stored in a matrix data: the first column is Y, 
    % the second column is X, and the third column is the category
    Y = data(:, 1);  
    X = data(:, 2);  
    categories = data(:, 3);  

    % Calculate Pearson correlation coefficient r and P-value
    [r, p] = corr(X, Y);

    % Perform linear regression fitting
    model = fitlm(X, Y);

    % Extract the slope and intercept of the fit
    coeffs = model.Coefficients.Estimate;  

    % The fitted linear equation
    x_fit = linspace(min(X), max(X), 100);  
    [y_fit, ci] = predict(model, x_fit');  

    % Specify category color
    custom_colors = [0, 0, 0; 67, 205, 128; 238, 99, 99] / 255;  

    % Draw a scatter plot
    fig = figure;
    fig.Position = [100, 100, 800, 600];

    % Set the entire graphic background to white or transparent
    %fig.Color = 'white';  
    fig.Color = 'none'; 

    hold on;
    unique_categories = unique(categories);
    for i = 1:length(unique_categories)
        % Filter the points of the current category
        current_category = unique_categories(i);
        idx = categories == current_category;
        
        % Draw scatter points for the current category
        scatter(X(idx), Y(idx), 40, 'MarkerFaceColor', custom_colors(i, :), ...
            'MarkerEdgeColor', 'k');
    end

    % Draw confidence intervals
    fill([x_fit, fliplr(x_fit)], [ci(:, 1)', fliplr(ci(:, 2)')], ...
        [211, 211, 211] / 255, 'EdgeColor', 'none', 'FaceAlpha', 0.5);  

    % Draw a fitted line
    plot(x_fit, y_fit, 'r-', 'LineWidth', 2);  

    hold off;


    set(gca, 'XTickLabel', [], 'YTickLabel', []);

    % Only retain the left Y-axis and bottom X-axis
    ax = gca;
    ax.LineWidth = 1.5;  
    ax.TickDir = 'in';   
    ax.XAxis.TickDirection = 'in';  
    ax.YAxis.TickDirection = 'in';  


    xticks(0:0.2:1);
    yticks(-1:0.2:1);
    xlim(x_lim);


    export_fig(out_name, '-m6', '-q400');


end




