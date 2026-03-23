clc;clear;

work_dir = '/home/jinlong/VDisk1/Jinlong/2025_HFR_ai_mod/';
ret_dir = '/home/jinlong/VDisk1/Jinlong/2025_HFR_ai_mod/result';

% network_label_all sleep_info sublist sub_stage_cell 
load([ret_dir '/IndiPar_net_multi_7/Network_ind_ret8min.mat']);
load([ work_dir '/sleep_subinfo.mat']);
load([ work_dir '/resource/fs_LR_32k_no_medial_wall_index.mat']);  %no_medial_wall_index
% line 6: eo1
sleep_sub_stage_cnt(:, ~filter_sub_idx) = [];
sleep_sub_file_info(:, ~filter_sub_idx) = [];
sleep_sub_stage_info(:, ~filter_sub_idx) = [];
sublist(~filter_sub_idx) = [];
network_label_all(:, ~filter_sub_idx) = [];

load([ret_dir '/IndiPar_net_multi_7/Network_ind_subcortex8min.mat']);
% addpath('/nd_disk3/guoyuan/Xinyu/plot_fig_subcortex/dependencies/export_fig/');

network_n = 7;
stageNames = {'Awake', 'N1', 'N2', 'N3', 'REM'};
rept_svm = 20;
kfold = 5;
nTrees = 100;
null_rept = 20;
balance_rpt_time = 5;

% create match table
% network_label_all = network_label_subcortex_lisa;
network_label_all = network_label_subcortex_lisa;
% feature_data = [];

template_dlabel = '/home/jinlong/VDisk1/Jinlong/2025_HFR_ai_mod/resource/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_netassignments_LR.dlabel.nii';

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
sequence_i = 1;
%%

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
                
                if mask_idx == 5 || mask_idx == 1
                    label_tmp(setdiff(1:64984, no_medial_wall_index)) = [];
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
    save(['/home/jinlong/VDisk1/Jinlong/2025_HFR_ai_mod/result/feature_data_' mask_name_list_w{mask_idx} '.mat'], 'feature_data');
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



% svm prediction 
sleep_stage_classify_acc = [];
sleep_stage_classify_confMatrix = cell(balance_rpt_time, 1);
sleep_stage_classify_beta = cell(balance_rpt_time, 1);
sleep_stage_classify_bias = cell(balance_rpt_time, 1);

allTrueLabels_list = [];
allPredictedScores_list = [];

externalConfusionMatrix_all = [];
nc_time_all = 0;

disp('process svm classification');

% if isempty(gcp('nocreate'))
%     parpool('local', 20);
% end

for balance_i = 1:balance_rpt_time
    disp(['process balance: ' num2str(balance_i) ]); tic;
    feature_data = feature_data_balance_cell{balance_i};
    
    [externalAccuracies, externalConfusionMatrix, beta_Matrix, bias_Matrix, allTrueLabels, allPredictedScores, nc_time] = svm_predict_stage(feature_data, rept_svm, kfold, ...
    '');
    sleep_stage_classify_acc = [sleep_stage_classify_acc; externalAccuracies];
    if isempty(externalConfusionMatrix_all)
        externalConfusionMatrix_all = externalConfusionMatrix;
    else
        externalConfusionMatrix_all = externalConfusionMatrix_all + externalConfusionMatrix;
    end
    nc_time_all = nc_time_all + nc_time;
    
    sleep_stage_classify_beta{balance_i} = beta_Matrix;
    sleep_stage_classify_bias{balance_i} = bias_Matrix;

    allTrueLabels_list = [allTrueLabels_list; allTrueLabels];
    allPredictedScores_list = [allPredictedScores_list; allPredictedScores];
%     draw_confMatrix(externalConfusionMatrix, stageNames, ...
%         ['/nd_disk3/guoyuan/Jinlong/sleep_ret/HFR_ret/pic/SVM_cbm_Mean_ConfMatrix_heatmap_Multiclass_network' num2str(network_n) '_balance' num2str(balance_i) '.png']);
		toc;
end

externalConfusionMatrix_all = externalConfusionMatrix_all ./ nc_time_all;
externalConfusionMatrix_all = externalConfusionMatrix_all ./ sum(externalConfusionMatrix_all, 2);

% ROC curve 
y_class_n = length(unique(allTrueLabels_list));
ROC_x = [];
ROC_y = [];
AUC_list = [];
for classIdx = 1:y_class_n
    
    [X, Y, ~, AUC] = perfcurve(allTrueLabels_list, allPredictedScores_list(:, classIdx), classIdx, 'XCrit', 'fpr', 'YCrit', 'tpr');
    
    ROC_x = [ROC_x; X];
    ROC_y = [ROC_y; Y];
    AUC_list = [AUC_list; AUC];
end

% [ROC_x_multi, ROC_y_multi, ~, AUC_multi] = perfcurve(allTrueLabels_list, allPredictedScores_list(sub2ind(size(allPredictedScores_list), (1:size(allPredictedScores_list, 1))', allTrueLabels_list)), 4, 'XCrit', 'fpr', 'YCrit', 'tpr');

X_macro = linspace(0, 1, 5);
Y_macro = zeros(size(X_macro));
for classIdx = 1:4
    X = ROC_x((1+(classIdx-1)*5):(classIdx*5));
    Y = ROC_y((1+(classIdx-1)*5):(classIdx*5));
    Y_interp = interp1(X, Y, X_macro, 'linear', 'extrap');  % ??D?2??¦Ì
    Y_macro = Y_macro + Y_interp;
end
Y_macro = Y_macro / 4;
ROC_x_multi = X_macro';
ROC_y_multi = Y_macro';
AUC_multi = mean(AUC_list);

% null model prediction

null_sleep_stage_classify_acc = [];



% random sleep stage label
for balance_i = 1:balance_rpt_time
	null_sleep_stage_classify_acc_tmp = zeros(null_rept, 1);
    new_feature_data = feature_data_balance_cell{balance_i};
    
    for rpti = 1:null_rept
        disp(['null-model svm predict. randi:' num2str(rpti)]); tic;
        feature_temp = new_feature_data;
        y = new_feature_data(:, size(new_feature_data, 2));
        feature_temp(:, size(new_feature_data, 2)) = y(randperm(length(y)));
        
        [externalAccuracies, externalConfusionMatrix, beta_Matrix, bias_Matrix, ~, ~, ~] = svm_predict_stage(feature_temp, 1, kfold, '');
        
        null_sleep_stage_classify_acc_tmp(rpti) = externalAccuracies;
        toc;
    end
    null_sleep_stage_classify_acc = [null_sleep_stage_classify_acc; null_sleep_stage_classify_acc_tmp];
end

% % % leave one network out prediction
% % leave1_sleep_stage_classify_acc = [];
% % leave1_sleep_stage_classify_acc_org = [];
% % leave1_sleep_stage_classify_confMatrix = cell(network_n, balance_rpt_time);
% % for balance_i = 1:balance_rpt_time
% %     leave1_sleep_stage_classify_acc_tmp = zeros(network_n, rept_svm);
% %     feature_data = feature_data_balance_cell{balance_i};
% %     
% %     for network_i = 1:network_n
% % %         if network_i == 4
% % %             continue;
% % %         end
% %         disp(['leave-one-out svm predict. exclude network:' num2str(network_i)]); tic;
% %         feature_temp = feature_data;
% %         feature_temp(:, network_i ) = [];
% %         [externalAccuracies, externalConfusionMatrix, beta_Matrix, bias_Matrix, allTrueLabels, allPredictedScores, nc_time] = svm_predict_stage(feature_data, rept_svm, kfold, ...
% %     '');
% %         leave1_sleep_stage_classify_acc_tmp(network_i, :) = externalAccuracies;
% %         leave1_sleep_stage_classify_confMatrix{network_i, balance_i} = externalConfusionMatrix;
% %         toc;
% %     end
% %     leave1_sleep_stage_classify_acc = [leave1_sleep_stage_classify_acc leave1_sleep_stage_classify_acc_tmp];
% % end

% delete(gcp('nocreate'));

save(['/home/jinlong/VDisk1/Jinlong/2025_HFR_ai_mod/out/SVM_onehot' num2str(sequence_i) '_all_multiclass_all_50_rpt_classify_network' num2str(network_n) '_' mask_name_list_w{mask_idx} '.mat' ], ...
 	'sleep_stage_classify_acc', 'externalConfusionMatrix_all', ...
     'sleep_stage_classify_beta', 'sleep_stage_classify_bias', ...
 	'null_sleep_stage_classify_acc', ...
 	'ROC_x',  'ROC_y', 'AUC_list', 'ROC_x_multi', 'ROC_y_multi', 'AUC_multi', 'allTrueLabels_list', 'allPredictedScores_list');


function [externalAccuracies, externalConfusionMatrix, beta_Matrix, bias_Matrix, allTrueLabels, allPredictedScores, nc_time] = svm_predict_stage(data, nRepeat, kfold_time, ROC_file_stem)
    % Assume X is the feature matrix, y is the label vector
    n = size(data, 2);
    X = data(:, 1:n-1); % Feature data
    y = data(:, n); % Label data
    y_class = unique(y);
    y_class_n = length(y_class);

    % Set external cross-validation parameters (2-fold cross-validation, repeated 100 times)
    k = kfold_time;  % Binary cross-validation (2-fold)

    % Grid search hyperparameter range (e.g., C and gamma parameters for SVM)
    BoxConstraints = [0.1, 1, 10];  % Candidate values for BoxConstraint
    KernelScales = [0.1, 1, 10];    % Candidate values for KernelScale

    % Store results of external cross-validation
    externalAccuracies = zeros(nRepeat, 1);
    externalConfusionMatrix = zeros(y_class_n, y_class_n);
    allTrueLabels = [];  % True labels
    allPredictedScores = [];
    allPredictedScores_multi = [];
    beta_Matrix = zeros(n-1, y_class_n, y_class_n);
    bias_Matrix = zeros(y_class_n, y_class_n);
    nc_time = 0;

    for repeat = 1:nRepeat
        % External cross-validation partition
        cv = cvpartition(length(y), 'KFold', k);
        
        % Store internal cross-validation accuracy (for hyperparameter tuning)
        foldAccuracies = zeros(cv.NumTestSets, 1);
        
        for fold = 1:cv.NumTestSets
            % Get indices for training and test sets
            contain_all_class_flg = 100;
            while contain_all_class_flg
                trainIdx = cv.training(fold);
                testIdx = cv.test(fold);
                
                % Get training and test set data
                XTrain = X(trainIdx, :);
                yTrain = y(trainIdx);
                XTest = X(testIdx, :);
                yTest = y(testIdx);
                if nnz(unique(yTrain) == y_class) == y_class_n && nnz(unique(yTest) == y_class) == y_class_n
                    break;
                end
                contain_all_class_flg = contain_all_class_flg - 1;
            end
            if contain_all_class_flg == 0
                break;
            end
            
            % One-vs-One strategy: train all possible binary SVMs
            classes = y_class;  % Get all class labels
            numClasses = length(classes);
            models = cell(numClasses, numClasses);  % Store SVM models for each pair of classes
            
            % Train SVM models for all class pairs
            for i = 1:numClasses-1
                for j = i+1:numClasses
                    % Get training data for the current class pair
                    idx = (yTrain == classes(i)) | (yTrain == classes(j));
                    XPair = XTrain(idx, :);
                    yPair = yTrain(idx);
                    
                    % Convert labels to binary labels (class i as 1, class j as -1)
                    yPair(yPair == classes(i)) = 1;
                    yPair(yPair == classes(j)) = -1;
                    
                    % Train SVM model
                    bestSVMModel = svm_model_optimal(XPair, yPair, BoxConstraints, KernelScales);
                    models{i,j} = bestSVMModel;
                    beta_Matrix(:, i, j) = beta_Matrix(:, i, j) + bestSVMModel.Beta(:);
                    bias_Matrix(i, j) = bias_Matrix(i, j) + bestSVMModel.Bias;
                end
            end
            
            % Make predictions on the test set
            yPred = zeros(length(yTest), 1);
            decisionValues = zeros(length(yTest), numClasses);  % Store decision values for each sample for each class
            for i = 1:length(yTest)
                % Calculate which class each sample belongs to
                votes = zeros(numClasses, 1);
                for iClass = 1:numClasses-1
                    for jClass = iClass+1:numClasses
                        % Get SVM model for the current class pair
                        model = models{iClass,jClass};
                        
                        % Predict current test sample
                        testSample = XTest(i, :);
                        prediction = predict(model, testSample);
                        % Update decision values (decision values for each class)
                        if prediction == 1
                            decisionValues(i, iClass) = decisionValues(i, iClass) + 1;
                        else
                            decisionValues(i, jClass) = decisionValues(i, jClass) + 1;
                        end
                        
                        % Increase votes based on SVM prediction result
                        if prediction == 1
                            votes(iClass) = votes(iClass) + 1;
                        else
                            votes(jClass) = votes(jClass) + 1;
                        end
                    end
                end
                
                % Select the class with the most votes
                [~, maxVoteClass] = max(votes);
                yPred(i) = classes(maxVoteClass);
            end

            % Add true labels and predicted decision values to storage list
            allTrueLabels = [allTrueLabels; yTest];
            allPredictedScores = [allPredictedScores; decisionValues];
            
            % Calculate accuracy for this cross-validation fold
            foldAccuracies(fold) = sum(yPred == yTest) / length(yTest);
            confMat = confusionmat(yTest, yPred);
            
            externalConfusionMatrix = externalConfusionMatrix + confMat;
            
            
            nc_time = nc_time + 1;
        end
        
        % Calculate external cross-validation accuracy
        externalAccuracies(repeat) = mean(foldAccuracies);
    end

    % Calculate final average accuracy and standard deviation
    avgAccuracy = mean(externalAccuracies);
    stdAccuracy = std(externalAccuracies);

    % Confusion matrix normalization (commented out)
%     externalConfusionMatrix = externalConfusionMatrix ./ nc_time;
%     externalConfusionMatrix = externalConfusionMatrix ./ sum(externalConfusionMatrix, 2);
%     externalConfusionMatrix = normalize(externalConfusionMatrix, 'norm', 1);
%     externalConfusionMatrix(isnan(externalConfusionMatrix)) = 0;

    beta_Matrix = beta_Matrix ./ nc_time;
    bias_Matrix = bias_Matrix ./ nc_time;

    fprintf('Average accuracy: %.4f\n', avgAccuracy);
    fprintf('Accuracy standard deviation: %.4f\n', stdAccuracy);

    allPredictedScores = allPredictedScores ./ 6;
    if ~isempty(ROC_file_stem)
        draw_ROC(allTrueLabels, allPredictedScores, {'Awake', 'N1', 'N2', 'N3', 'REM'}, ROC_file_stem);
    end

end

function bestSVMModel = svm_model_optimal(XTrain, yTrain, BoxConstraints, KernelScales)
    % Perform grid search for hyperparameter tuning
    bestAccuracy = -Inf;
    bestBoxConstraint = NaN;
    bestKernelScale = NaN;
    
    % Iterate through each hyperparameter combination
    for boxC = BoxConstraints
        for kernelS = KernelScales
            % Train SVM model using current BoxConstraint and KernelScale
            SVMModel = fitcsvm(XTrain, yTrain, 'KernelFunction', 'linear', ...
                'BoxConstraint', boxC, 'KernelScale', kernelS, 'CrossVal', 'on');
            
            % Get cross-validation accuracy for the current model
            cvAccuracy = 1 - kfoldLoss(SVMModel);  % Calculate accuracy
            
            % Record optimal hyperparameters
            if cvAccuracy > bestAccuracy
                bestAccuracy = cvAccuracy;
                bestBoxConstraint = boxC;
                bestKernelScale = kernelS;
            end
        end
    end
    
    % Train the final model using optimal hyperparameters
    bestSVMModel = fitcsvm(XTrain, yTrain, 'KernelFunction', 'linear', ...
        'BoxConstraint', bestBoxConstraint, 'KernelScale', bestKernelScale);
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

	% Define colors: gradient from dark gray to red to dark red
    colorMap = [
        0.2, 0.2, 0.2;    % Dark gray
        0.5, 0.5, 0.5;    % Gray
        0.9, 0.9, 0.8;    % Beige
        1, 0.7, 0.7;      % Pink
        1, 0, 0;          % Red
        0.6, 0, 0;        % Dark red
    ];

	% Refine color map to achieve gradient effect
    nSteps = 256;  % Increase detail, generate more color transitions
    gradients = linspace(0, 1, nSteps);  % Create gradient values from 0 to 1

    % Interpolate color map
    smoothColorMap = [];
    for i = 1:size(colorMap, 1) - 1
        % Get RGB values of each adjacent color pair
        startColor = colorMap(i, :);
        endColor = colorMap(i + 1, :);
        
        % Interpolate between these two colors to generate gradient colors
        gradColors = [linspace(startColor(1), endColor(1), nSteps/size(colorMap, 1))', ...
                    linspace(startColor(2), endColor(2), nSteps/size(colorMap, 1))', ...
                    linspace(startColor(3), endColor(3), nSteps/size(colorMap, 1))'];
        
        % Add gradient colors to smooth color map
        smoothColorMap = [smoothColorMap; gradColors];
    end

    % Set color map
    colormap(smoothColorMap);  % Set gradient color map

	caxis([0 0.6]);     % Set color range
	title('Confusion matrix');
	xlabel('Predicted', 'FontSize', 18);
	ylabel('Ground truth', 'FontSize', 18);

	% Set x-axis and y-axis ticks and labels
	xticks(1:4);                     % Set x-axis ticks from 1 to 4
	yticks(1:4);                     % Set y-axis ticks from 1 to 4
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
	caxis([0 0.6]);  % Set colorbar value range same as caxis

	% Set colorbar ticks (display values), step size 0.2
	c.Ticks = [0, 0.25, 0.5, 0.85, 1];  % Set colorbar ticks
	c.TickLabels = {'0', '0.3', '0.6', '0.85', '1'};  % Display corresponding values

	% Set colorbar position
	c.Position = [0.91, 0.2, 0.03, 0.6];  % Adjust position to avoid overlap with plot

% 	saveas(gcf, out_name);
    export_fig(out_name, '-m6', '-q100');
end


function draw_ROC(allTrueLabels, allPredictedScores, stage_Names, ROC_file_stem)
    % Plot ROC curves
    fig = figure;
    hold on;
    y_class_n = length(unique(allTrueLabels));
    for classIdx = 1:y_class_n
        % Calculate ROC curve for each class
        [X, Y, ~, AUC] = perfcurve(allTrueLabels, allPredictedScores(:, classIdx), classIdx, 'XCrit', 'fpr', 'YCrit', 'tpr');
        plot(X, Y, 'LineWidth', 2, 'DisplayName', [stage_Names{classIdx} ' (AUC = ' num2str(AUC) ')']);
    end
    xlabel('False Positive Rate');
    ylabel('True Positive Rate');
    title('ROC Curves for Multiclass');
    legend('Location', 'southeast');  % Place legend in bottom right corner
    grid on;
    hold off;
%     saveas(gcf, [ROC_file_stem '.png']);
    export_fig([ROC_file_stem '.png'], '-m6', '-q100');
    close(fig);
end

function onehot_mat = onehot_encode(label, n_class)

    n_sample = length(label);
    
    onehot_mat = full(sparse(1:n_sample, label, 1, n_sample, n_class));
end