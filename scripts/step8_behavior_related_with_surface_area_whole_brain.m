% Author: Jinlong Li and Guoyuan Yang, BIT.
% Fit the correlation curve between the percentage of network size 
% and PSQI for a single network
% 
clc;clear;

Ancillary_data_file = '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/resource/Ancillary_data_final.xlsx';
Model_sleep_data_file = '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/resource/Model_sleep_wake_SWA_N2562.txt';

% network_label_all sleep_info sublist sub_stage_cell 
load('/nd_disk3/guoyuan/Jinlong/sleep_ret/HFR_ret/IndiPar_net7/Network_ind_ret8min.mat');
load('/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/resource/fs_LR_32k_no_medial_wall_index.mat');  %no_medial_wall_index
info = readtable(Ancillary_data_file);
work_dir = '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/';
load([ work_dir '/sleep_subinfo.mat']);
load( '/nd_disk3/guoyuan/Jinlong/sleep_ret/HFR_ret/IndiPar_net7/Network_ind_subcortex_cleanup8min.mat');


% line 6: eo1

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
mask_idx = 1;

network_n = 7;
vertex_n = 64984;

network_size_ret = cell(size(network_label_subcortex_lisa));



for stai = 1:5
    
    nc = 0;
    for subi = 1:size(network_label_subcortex_lisa, 2)
        if isempty(network_label_subcortex_lisa{stai, subi})
            continue;
        end
        temp_size_sta = zeros(1, network_n);

            label_tmp = network_label_subcortex_lisa{stai, subi};
            label_tmp = label_tmp(mask_list_w{mask_idx});
            for ni = 1:network_n
                idx = find(label_tmp == ni);
                if isempty(idx)
                    disp([num2str(stai) ' ' num2str(subi)])
                end
                temp_size_sta(ni) = length(idx) ;
            end

        network_size_ret{stai, subi} = temp_size_sta ./ length(label_tmp);
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
    for st = 1:4
        if isempty(network_size_ret{st+1, subi})
            network_size_ret1{st, subi} = [];
            continue;
        end
        network_size_ret1{st, subi} = network_size_ret{st+1, subi} - network_size_ret{1, subi};
    end
end


sub1 = info.ID;
sub2 = zeros(1, length(sub1));
for ki = 1:length(sub1)
    ta = strsplit(sub1{ki}, 'sub'); sub2(ki) = str2double(ta{2});
end
[ia, ib] = ismember(sublist, sub2);
sub3 = sublist(ia);
sub4 = sub2(ib(ia));
info = info(ib(ia), :);



behavior_list = {info.PSQI, info.ESS, info.MEQ, info.FSS};
bahevior_name_list = {'PSQI', 'ESS', 'MEQ', 'FSS'};
if network_n == 7
    label_name = {'Visual','Somatomotor','Dorsal attention', 'Action mode', 'Limbic', 'Frontoparietal', 'Default'};
elseif network_n == 17
    label_name = {'Auditory','Dorsal Attention A','Control A', 'Somatomotor A', 'Salience/VenAttn B', ...
'Default B', 'Default C', 'Visual C', 'Visual A', 'Dorsal Attention B', 'Temoral Pariental', 'Control B','Visual B',...
'Control C','Default A','Salience/VenAttn A', 'Somatomotor B'};
end
stageNames = {'Awake', 'N1', 'N2', 'N3', 'REM'};

network_size_variability = cell(network_n, length(behavior_list));

P_FDR_list = cell(length(behavior_list), 1);
P_org_list = cell(length(behavior_list), 1);


for bei = 1:1%length(behavior_list)
    p_list_org = zeros(network_n, 3);
    for networki = 1:network_n
        stage_size_variability = cell(3, 1);
        for st = 1:4
            kidx = 0;
            temp_size_variability = zeros(10, 2);
            
            for subi = 1:length(sub3)
                if isempty(network_size_ret1{st, subi})
                    continue;
                end
                kidx = kidx + 1;
                temp_size_variability(kidx, :) = [network_size_ret1{st, subi}(networki) behavior_list{bei}(subi)];
            end

            stage_size_variability{st} = temp_size_variability;
            [r, p] = corr(temp_size_variability(:, 2), temp_size_variability(:, 1));
            p_list_org(networki, st) = p;
  
        end
        network_size_variability{networki, bei} = stage_size_variability;
    end
    P_org_list{bei} = p_list_org;
    p_list_org = reshape(p_list_org, [], 1);
    [pthr,pcor,padj] = fdr(p_list_org);
    p_list_mdr = reshape(padj, network_n, 4);
    P_FDR_list{bei} = p_list_mdr;

end

x_tick_cell = {
  0:4:12, 0:4:12, 0:5:10, 0:4:12
};
y_tick_cell = {
  -0.1:0.1:0.1, -0.14:0.04:0.12, -0.1:0.05:0.1, -0.12:0.04:0.12  
};
% -0.1 0.1
for bei = 1:1%length(behavior_list)
    p_list_fdr = P_FDR_list{bei};
    for networki = [ 2 4]
        stage_size_variability = network_size_variability{networki, bei};
        for st = 3:3
            temp_size_variability = stage_size_variability{st};
            draw_plot(temp_size_variability, bahevior_name_list{bei}, label_name{networki}, stageNames{st+1}, ['/nd_disk3/guoyuan/Jinlong/sleep_ret/HFR_ret/pic/behavior_related_' num2str(network_n)  'net_' bahevior_name_list{bei} '_net' num2str(networki) '_st' num2str(st)  '.png'], P_org_list{bei}(networki, st), x_tick_cell{st}, y_tick_cell{st});
            draw_plot_clean(temp_size_variability, bahevior_name_list{bei}, label_name{networki}, stageNames{st+1}, ['/nd_disk3/guoyuan/Jinlong/sleep_ret/HFR_ret/pic/behavior_related_' num2str(network_n)  'net_' bahevior_name_list{bei} '_net' num2str(networki) '_st' num2str(st)  '_clean.png'], P_org_list{bei}(networki, st), x_tick_cell{st}, y_tick_cell{st});
        end
        
    end

end



    


function draw_plot(data, behavior_name, network_name, stage_name, out_name, p_FDR, x_tick_list, y_tick_list)
    
    % Assuming the data is stored in a matrix data: 
    % the first column is Y, and the second column is X


    Y = data(:, 1);  
    X = data(:, 2);  

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
    custom_colors = [0, 134, 139; 67, 205, 128; 238, 99, 99] / 255;  

    % Draw a scatter plot
    fig = figure;
    fig.Position = [100, 100, 800, 600];

    % Set the entire graphic background to white or transparent
    %fig.Color = 'white';  
    fig.Color = 'none'; 

    scatter(X, Y, 40, 'MarkerFaceColor', 'k', ...
            'MarkerEdgeColor', 'k');
    hold on;


    % Draw confidence intervals
    fill([x_fit, fliplr(x_fit)], [ci(:, 1)', fliplr(ci(:, 2)')], ...
        [211, 211, 211] / 255, 'EdgeColor', 'none', 'FaceAlpha', 0.5);  

    % Draw a fitted line
    plot(x_fit, y_fit, 'r-', 'LineWidth', 2);  
    hold off;

    

    % Add labels
    xlabel(behavior_name, 'FontSize', 16);
    ylabel('variation of surface area', 'FontSize', 16);
    title(['Linear fit: ' network_name '(' stage_name ' vs Awake)'], 'FontSize', 16);

    set(gca, 'FontSize', 16);       
    set(gca, 'XColor', 'k', 'YColor', 'k', 'Box', 'off');  
    
    % Only retain the left Y-axis and bottom X-axis
    ax = gca;
    ax.LineWidth = 1.5;  
    ax.TickDir = 'in';   
    ax.XAxis.TickDirection = 'in';  
    ax.YAxis.TickDirection = 'in';  

    % Mark the values of r and P in the bottom right corner

    text(max(X), min(Y), sprintf('r = %.2f\nP = %.3g', r, p_FDR), 'FontSize', 16, 'Color', 'black', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

    xticks(x_tick_list);
    yticks(y_tick_list);

    export_fig(out_name, '-m6', '-q400');
    close(fig)
end




function draw_plot_clean(data, behavior_name, network_name, stage_name, out_name, p_FDR, x_tick_list, y_tick_list)
    
    % Assuming the data is stored in a matrix data: 
    % the first column is Y, and the second column is X


    Y = data(:, 1);  
    X = data(:, 2);  

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
    custom_colors = [0, 134, 139; 67, 205, 128; 238, 99, 99] / 255;  

    % Draw a scatter plot
    fig = figure;
    fig.Position = [100, 100, 800, 600];

    % Set the entire graphic background to white or transparent
    %fig.Color = 'white';  
    fig.Color = 'none'; 

    scatter(X, Y, 40, 'MarkerFaceColor', 'k', ...
            'MarkerEdgeColor', 'k');
    hold on;


    % Draw confidence intervals
    fill([x_fit, fliplr(x_fit)], [ci(:, 1)', fliplr(ci(:, 2)')], ...
        [211, 211, 211] / 255, 'EdgeColor', 'none', 'FaceAlpha', 0.5);  

    % Draw a fitted line
    plot(x_fit, y_fit, 'r-', 'LineWidth', 2);  
    hold off;

    

    % Add labels

    set(gca, 'FontSize', 16);       
    set(gca, 'XColor', 'k', 'YColor', 'k', 'Box', 'off');  
    set(gca, 'XTickLabel', [], 'YTickLabel', []);

    ax = gca;
    ax.LineWidth = 1.5;  
    ax.TickDir = 'in';   
    ax.XAxis.TickDirection = 'in';  
    ax.YAxis.TickDirection = 'in';  


    xticks(x_tick_list);
    yticks(y_tick_list);

    export_fig(out_name, '-m6', '-q400');
    close(fig)
end




