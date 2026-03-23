clc;clear;

work_dir = '/home/jinlong/VDisk1/Jinlong/2025_HFR_ai_mod/';
ret_dir = '/home/jinlong/VDisk1/Jinlong/2025_HFR_ai_mod/result';
out_dir = '/home/jinlong/VDisk1/Jinlong/2025_HFR_ai_mod/out';

% network_label_all sleep_info sublist sub_stage_cell 
load([ret_dir '/IndiPar_net_multi_7/Network_ind_subcortex_cleanup8min.mat']);
load([ work_dir '/sleep_subinfo.mat']);
load([ work_dir '/resource/fs_LR_32k_no_medial_wall_index.mat']);  %no_medial_wall_index



% line 6: eo1
sleep_sub_stage_cnt(:, ~filter_sub_idx) = [];
sleep_sub_file_info(:, ~filter_sub_idx) = [];
sleep_sub_stage_info(:, ~filter_sub_idx) = [];
sublist(~filter_sub_idx) = [];
%network_label_all(:, ~filter_sub_idx) = [];

network_n = 7;
feature_label_n = 59412*7;
stageNames = {'Awake', 'N1', 'N2', 'N3', 'REM'};
rept_svm = 20;
kfold = 5;
balance_rpt_time = 5;


% roi mask
template_dlabel = [work_dir '/resource/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_netassignments_LR.dlabel.nii'];

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
mask_list_w = {cortex_mask_lh | cortex_mask_rh, cerebellum_mask_lh | cerebellum_mask_rh, thalamus_mask_lh | thalamus_mask_rh, striatum_mask_lh | striatum_mask_rh, ...
    cortex_mask_lh | cortex_mask_rh | cerebellum_mask_lh | cerebellum_mask_rh | thalamus_mask_lh | thalamus_mask_rh | striatum_mask_lh | striatum_mask_rh};
mask_name_list_w = {'cortex', 'cerebellum', 'thalamus', 'striatum', 'wholebrain'};


%% notify: 
% Modify the mask_idx variable to run results for different brain regions separately.
mask_idx = 1;

% Using network topology as features for classification results in a large number of features, 
% which leads to high computational time. 
% Parallel processing was implemented by distributing the workload across different subroutines 
% using sequence_i.If not needed, it can be ignored.
sequence_i = 3;
rng(sequence_i + 5904, 'twister');
%%

% create match table
network_label_all = network_label_subcortex_lisa;



nc_time1 = 0;

exist_flg == 0

if exist_flg == 0
    for stagei = 1:4
        for subi = 1:size(network_label_all, 2)
            if isempty(network_label_all{stagei, subi})
                continue;
            end
            for sessi = 1:length(network_label_all{stagei, subi})
                tic;
                label_tmp = network_label_all{stagei, subi}{sessi};
                label_tmp(mask_list_w{mask_idx} == 0) = [];
                if mask_idx == 1 || mask_idx == 5
                    label_tmp = label_tmp(no_medial_wall_index);
                end
                network_label_list_onehot = onehot_encode(label_tmp, 7);
                network_label_list_onehot = reshape(network_label_list_onehot, 1, []);
                
                network_label_list_onehot(end + 1) = stagei;
                nc_time1 = nc_time1 + 1;
                if ~exist('feature_data', 'var')
                    feature_data = network_label_list_onehot;
                else
                    feature_data(nc_time1, :) =  network_label_list_onehot;
                end
                toc;
            end
        end
    end
    save([out_dir '/SVM_onehot_binary_featuredata_mask' num2str(mask_idx) '.mat'], 'feature_data');
else
    load(['/home/jinlong/VDisk1/Jinlong/2025_HFR_ai_mod/result/feature_data_' mask_name_list_w{mask_idx} '.mat']);
end

% balance feature
num_ths = 178;
feature_data_balance_cell = cell(balance_rpt_time, 1);
for balance_i = 1:balance_rpt_time
    feature_data_balance = [];
    for stagei = 1:4
        nidx = find(feature_data(:, end) == stagei);
        r_nidx = randperm(length(nidx), num_ths);
        nidx_f = nidx(r_nidx);
        feature_data_balance = [feature_data_balance; feature_data(nidx_f, :)];
    end
    feature_data_balance_cell{balance_i} = feature_data_balance;
end

feature_label_n = size(feature_data_balance, 2) - 1;

% svm prediction 
sleep_stage_classify_acc = cell(5, 5);
sleep_stage_classify_std = zeros(5, 5);
sleep_stage_classify_confMatrix = cell(5, 5);
sleep_stage_classify_beta = zeros(5, 5, feature_label_n);
sleep_stage_classify_bias = zeros(5, 5);
for stagei = 1:3
	for stagej = stagei+1:4
        acc_list_tmp = [];
        for balance_i = sequence_i:sequence_i%balance_rpt_time
            feature_data = feature_data_balance_cell{balance_i};
            disp(['svm predict - ' num2str(balance_i) ' : ' stageNames{stagei} ' vs ' stageNames{stagej}]); tic;
		    feature_temp = feature_data(:, 1:feature_label_n);
		    label_bi = feature_data(:, feature_label_n+1) == stagei;
    
		    % filter sleep stage
		    label_bj = feature_data(:, feature_label_n+1) == stagej;
		    label_filter = (label_bi + label_bj) == 0;
		    feature_temp = [feature_temp label_bi];
		    feature_temp(label_filter, :) = [];
    
		    [externalAccuracies, externalConfusionMatrix, beta_Matrix, bias_Matrix] = svm_predict_stage(feature_temp, rept_svm, kfold);
		    acc_list_tmp = [acc_list_tmp; externalAccuracies];

		%draw_confMatrix(externalConfusionMatrix, {stageNames{stagei}, stageNames{stagej}}, ...
		%	['/nd_disk3/guoyuan/Jinlong/sleep_ret/HFR_ret/pic/Mean_ConfMatrix_heatmap_network' num2str(network_n) '_' stageNames{stagei} '_vs_' stageNames{stagej} '.png']);
		    toc;
        end
        sleep_stage_classify_acc{stagei, stagej} = acc_list_tmp;
% 	    sleep_stage_classify_confMatrix{stagei, stagej} = externalConfusionMatrix;
	    sleep_stage_classify_beta(stagei, stagej, :) = beta_Matrix(:);
% 	    sleep_stage_classify_bias(stagei, stagej) = bias_Matrix;
	end
end

% null model prediction
null_rept = 40;
null_sleep_stage_classify_acc = zeros(5, 5, null_rept);
null_sleep_stage_classify_std = zeros(5, 5, null_rept);
null_sleep_stage_classify_confMatrix = cell(5, 5, null_rept);
for stagei = 1:3
	for stagej = stagei+1:4
        for balance_i = sequence_i:sequence_i%balance_rpt_time
            tic;
            feature_data = feature_data_balance_cell{balance_i};
		    feature_temp = feature_data(:, 1:(size(feature_data, 2)-1));
		    label_bi = feature_data(:, size(feature_data, 2)) == stagei;
    
		    % filter sleep stage
		    label_bj = feature_data(:, size(feature_data, 2)) == stagej;
		    label_filter = (label_bi + label_bj) == 0;
		    feature_temp = [feature_temp label_bi];
		    feature_temp(label_filter, :) = [];
		    
		    % random sleep stage label
		    for rpti = 1:null_rept
			    disp(['null-model svm predict - ' num2str(balance_i) ': ' stageNames{stagei} ' vs ' stageNames{stagej} ' randi:' num2str(rpti)]); 
			    y = feature_temp(:, (size(feature_data, 2)));
			    feature_temp(:, size(feature_data, 2)) = y(randperm(length(y)));
			    
			    [externalAccuracies, externalConfusionMatrix, beta_Matrix, bias_Matrix] = svm_predict_stage(feature_temp, 1, kfold);
			    
			    null_sleep_stage_classify_acc(stagei, stagej, rpti+(balance_i-1)*null_rept) = externalAccuracies;
			    
            end
            toc;
        end
	end
end

% leave one network out prediction
leave1_sleep_stage_classify_acc = zeros(5, 5, network_n);
leave1_sleep_stage_classify_std = zeros(5, 5, network_n);
leave1_sleep_stage_classify_confMatrix = cell(5, 5, network_n);


save([out_dir '/SVM_onehot' num2str(sequence_i) '_40rpt_binary_' mask_name_list_w{mask_idx} '_classify_network' num2str(network_n) '.mat' ], ...
	'sleep_stage_classify_acc', 'sleep_stage_classify_confMatrix', ...
    'sleep_stage_classify_beta', 'sleep_stage_classify_bias', ...
	'null_sleep_stage_classify_acc', 'null_sleep_stage_classify_std', 'null_sleep_stage_classify_confMatrix');


function [externalAccuracies, externalConfusionMatrix, beta_Matrix, bias_Matrix] = svm_predict_stage(data, nRepeat, kfold_time)
    % Assume X is the feature matrix, y is the label vector
    n = size(data, 2);
    X = data(:, 1:n-1); % Your feature data
    y = data(:, n); % Your label data

    % Set external cross-validation parameters (2-fold cross-validation, repeated 100 times)
    % nRepeat = 100;  % Number of cross-validation repetitions
    k = kfold_time;  % Binary cross-validation (2-fold)

    % Grid search hyperparameter range (e.g., C and gamma parameters for SVM)
    BoxConstraints = [0.1, 1, 10];  % Candidate values for BoxConstraint
    KernelScales = [0.1, 1, 10];    % Candidate values for KernelScale

    % Store results of external cross-validation
    externalAccuracies = zeros(nRepeat, 1);
    externalConfusionMatrix = zeros(2, 2);
    beta_Matrix = zeros(n-1, 1);
    bias_Matrix = 0;
    nc_time = 0;

    for repeat = 1:nRepeat
        % External cross-validation partition
        cv = cvpartition(length(y), 'KFold', k);
        
        % Store internal cross-validation accuracy (for hyperparameter tuning)
        foldAccuracies = zeros(cv.NumTestSets, 1);
        
        for fold = 1:cv.NumTestSets
            % Get indices for training and test sets
            trainIdx = cv.training(fold);
            testIdx = cv.test(fold);
            
            % Get training and test set data
            XTrain = X(trainIdx, :);
            yTrain = y(trainIdx);
            XTest = X(testIdx, :);
            yTest = y(testIdx);
            
            % Perform grid search for hyperparameter tuning
            bestAccuracy = -Inf;
            bestBoxConstraint = NaN;
            bestKernelScale = NaN;
            
            % Iterate through each hyperparameter combination
% %           for boxC = BoxConstraints
% %               for kernelS = KernelScales
% %                   % Train SVM model using current BoxConstraint and KernelScale
% %                   SVMModel = fitcsvm(XTrain, yTrain, 'KernelFunction', 'linear', ...
% %                       'BoxConstraint', boxC, 'KernelScale', kernelS, 'CrossVal', 'on');
% %                   
% %                   % Get cross-validation accuracy for the current model
% %                   cvAccuracy = 1 - kfoldLoss(SVMModel);  % Calculate accuracy
% %                   
% %                   % Record optimal hyperparameters
% %                   if cvAccuracy > bestAccuracy
% %                       bestAccuracy = cvAccuracy;
% %                       bestBoxConstraint = boxC;
% %                       bestKernelScale = kernelS;
% %                   end
% %               end
% %           end
            
            % Train the final model using optimal hyperparameters
            bestSVMModel = fitcsvm(XTrain, yTrain, 'KernelFunction', 'linear');
            
            % Evaluate the model on the test set
            yPred = predict(bestSVMModel, XTest);
            foldAccuracies(fold) = sum(yPred == yTest) / length(yTest);
            confMat = confusionmat(yTest, yPred);
            
            externalConfusionMatrix = externalConfusionMatrix + confMat;
            beta_Matrix = beta_Matrix + bestSVMModel.Beta;
            bias_Matrix = bias_Matrix + bestSVMModel.Bias;

            nc_time = nc_time + 1;
        end
        
        % Calculate external cross-validation accuracy
        externalAccuracies(repeat) = mean(foldAccuracies);
    end

    % Calculate final average accuracy and standard deviation
    avgAccuracy = mean(externalAccuracies);
    stdAccuracy = std(externalAccuracies);

    % Confusion matrix normalization
    externalConfusionMatrix = externalConfusionMatrix ./ nc_time;
    externalConfusionMatrix = externalConfusionMatrix ./ sum(externalConfusionMatrix, 2);

    beta_Matrix = beta_Matrix ./ nc_time;
    bias_Matrix = bias_Matrix ./ nc_time;

    fprintf('Average accuracy: %.4f\n', avgAccuracy);
    fprintf('Accuracy standard deviation: %.4f\n', stdAccuracy);

end

function draw_confMatrix(confMatrix, stageNames, out_name)

    pMatrix = confMatrix;
    % Set significance levels
    alpha1 = 0.001;
    alpha2 = 0.01;
    alpha3 = 0.05;

    % Create color intensity plot
    fig = figure;
    fig.Position = [100, 100, 1000, 800];
    imagesc(pMatrix); % Visualize using imagesc

    % Define colors, from gray to beige to red, divided into 5 color blocks (one interval per 0.2)
    colorMap = [
        0.5, 0.5, 0.5;    % 0-0.2: Gray
        0.7, 0.7, 0.65;   % 0.2-0.4: Light gray to beige
        0.9, 0.9, 0.8;    % 0.4-0.6: Beige
        1, 0.7, 0.7;      % 0.6-0.8: Beige to red
        1, 0, 0;          % 0.8-1.0: Red
    ];

    % Map to values from 0 to 1
    colormap(colorMap);  % Set color map

    caxis([0 1]);     % Set color range
    title('Confusion matrix');
    xlabel('Predicted', 'FontSize', 18);
    ylabel('Ground truth', 'FontSize', 18);

    % Set x-axis and y-axis ticks and labels
    xticks(1:2);                     % Set x-axis ticks from 1 to 2
    yticks(1:2);                     % Set y-axis ticks from 1 to 2
    xticklabels(stageNames);         % Set x-axis tick labels
    yticklabels(stageNames);         % Set y-axis tick labels

    set(gca, 'FontSize', 18);       % Set axis tick font size
    set(gca, 'xticklabel', stageNames, 'FontSize', 18); % Set x-axis tick label font size
    set(gca, 'yticklabel', stageNames, 'FontSize', 18); % Set y-axis tick label font size

    % Add grid and significance value annotations
    hold on;
    for i = 1:size(pMatrix, 1)
        for j = 1:size(pMatrix, 2)
            % Display values, unified black font
            text(j, i, sprintf('%.3f', pMatrix(i, j)), ...
                'HorizontalAlignment', 'center', ...
                'Color', 'black', 'FontSize', 12, 'FontWeight', 'bold');
        end
    end
    % Manually draw black borders for each grid cell
    for i = 1:size(pMatrix, 1)
        for j = 1:size(pMatrix, 2)
            rectangle('Position', [j-0.5, i-0.5, 1, 1], ...
                    'EdgeColor', 'k', 'LineWidth', 1); % Draw black border
        end
    end
    hold off;

    % Add colorbar
    c = colorbar;  % Display colorbar

    % Remove colorbar ticks
    c.Ticks = [];  % Hide ticks
    c.TickLabels = {};  % Hide tick labels

    % Set colorbar ticks (display values)
    caxis([0 1]);  % Set colorbar value range same as caxis

    % Set colorbar ticks (display values), step size 0.2
    c.Ticks = 0:0.2:1;  % Set colorbar ticks, one every 0.2
    c.TickLabels = arrayfun(@(x) sprintf('%.1f', x), 0:0.2:1, 'UniformOutput', false);  % Display corresponding values

    % Set colorbar position
    c.Position = [0.91, 0.2, 0.03, 0.6];  % Adjust position to avoid overlap with plot

    saveas(gcf, out_name);
end

function onehot_mat = onehot_encode(label, n_class)

    n_sample = length(label);
    
    onehot_mat = full(sparse(1:n_sample, label, 1, n_sample, n_class));
end