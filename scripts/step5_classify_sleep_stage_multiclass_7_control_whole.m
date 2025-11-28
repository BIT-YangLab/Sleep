% Author: Jinlong Li and Guoyuan Yang, BIT.
% multiple classification for sleep stage with network size ratio by SVM
% construct multiple classify model of SVM by one-against-one measure

clc;clear;

work_dir = '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/';
ret_dir = '/nd_disk3/guoyuan/Jinlong/sleep_ret/HFR_ret/';

% network_label_all sleep_info sublist sub_stage_cell 
load([ret_dir '/IndiPar_net_multi_7/Network_ind_ret5min.mat']);
load([ work_dir '/sleep_subinfo.mat']);
load([ work_dir '/resource/fs_LR_32k_no_medial_wall_index.mat']);  %no_medial_wall_index
% line 6: eo1
sleep_sub_stage_cnt(:, ~filter_sub_idx) = [];
sleep_sub_file_info(:, ~filter_sub_idx) = [];
sleep_sub_stage_info(:, ~filter_sub_idx) = [];
sublist(~filter_sub_idx) = [];
network_label_all(:, ~filter_sub_idx) = [];

load([ret_dir '/IndiPar_net_multi_7/Network_ind_subcortex_cleanup5min.mat']);
addpath('/nd_disk3/guoyuan/Xinyu/plot_fig_subcortex/dependencies/export_fig/');

network_n = 7;
stageNames = {'Awake', 'N1', 'N2', 'N3', 'REM'};
rept_svm = 40;
kfold = 5;
nTrees = 100;
null_rept = 50;
balance_rpt_time = 5;

% create match table
% network_label_all = network_label_subcortex_lisa;
network_label_all = network_label_subcortex_lisa;


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
mask_list_w = {cortex_mask_lh | cortex_mask_rh, cerebellum_mask_lh | cerebellum_mask_rh, thalamus_mask_lh | thalamus_mask_rh, striatum_mask_lh | striatum_mask_rh, ...
    cortex_mask_lh | cortex_mask_rh | cerebellum_mask_lh | cerebellum_mask_rh | thalamus_mask_lh | thalamus_mask_rh | striatum_mask_lh | striatum_mask_rh};
mask_name_list_w = {'cortex', 'cerebellum', 'thalamus', 'striatum', 'wholebrain'};

%% select different brain region
mask_idx = 5;

% control_sub_list_select control_sub_list
load([work_dir '/subject_control_select_7.mat']);

for random_i = 1:size(control_sub_list, 2)
    feature_data = [];
    for stagei = 1:4
	    for subi = 1:size(network_label_all, 2)
		    if isempty(network_label_all{stagei, subi})
			    continue;
		    end
            if nnz(control_sub_list_select{random_i} == subi) == 0
                continue;
            end
		    for sessi = 1:1%length(network_label_all{stagei, subi})
			    label_tmp = network_label_all{stagei, subi}{sessi};
                label_tmp = label_tmp(mask_list_w{mask_idx});
			    network_size_list = zeros(1, network_n+1);
			    for networki = 1:network_n
				    idx = find(label_tmp == networki);
				    network_size_list(networki) = length(idx) / nnz(mask_list_w{mask_idx});
			    end
			    network_size_list(networki + 1) = stagei;
			    feature_data = [feature_data; network_size_list];
		    end
	    end
    end
    feature_data_balance_cell{random_i} = feature_data;
end



% svm prediction 
sleep_stage_classify_acc = [];
sleep_stage_classify_confMatrix = cell(balance_rpt_time, 1);
sleep_stage_classify_beta = cell(balance_rpt_time, 1);
sleep_stage_classify_bias = cell(balance_rpt_time, 1);

allTrueLabels_list = [];
allPredictedScores_list = [];

nc_time_all = 0;

for balance_i = 1:balance_rpt_time
    feature_data = feature_data_balance_cell{balance_i};
    [externalAccuracies, externalConfusionMatrix, beta_Matrix, bias_Matrix,  allTrueLabels, allPredictedScores, nc_time] = svm_predict_stage(feature_data, rept_svm, kfold, ...
    '');
    sleep_stage_classify_acc = [sleep_stage_classify_acc; externalAccuracies];
    sleep_stage_classify_confMatrix{balance_i} = externalConfusionMatrix;
    sleep_stage_classify_beta{balance_i} = beta_Matrix;
    sleep_stage_classify_bias{balance_i} = bias_Matrix;

    allTrueLabels_list = [allTrueLabels_list; allTrueLabels];
    allPredictedScores_list = [allPredictedScores_list; allPredictedScores];

    if balance_i == 1
        externalConfusionMatrix_list = externalConfusionMatrix;
    else
        externalConfusionMatrix_list = externalConfusionMatrix_list + externalConfusionMatrix;
    end
    nc_time_all = nc_time_all + nc_time;


end

externalConfusionMatrix_list = externalConfusionMatrix_list ./ nc_time_all;

% ROC curve 
y_class_n = length(unique(allTrueLabels_list));
ROC_x = [];
ROC_y = [];
AUC_list = [];
for classIdx = 1:y_class_n
    % Calculate the ROC curve for each category
    [X, Y, ~, AUC] = perfcurve(allTrueLabels_list, allPredictedScores_list(:, classIdx), classIdx, 'XCrit', 'fpr', 'YCrit', 'tpr');
    
    ROC_x = [ROC_x; X];
    ROC_y = [ROC_y; Y];
    AUC_list = [AUC_list; AUC];
end

% null model prediction

null_sleep_stage_classify_acc = [];



% random sleep stage label
for balance_i = 1:balance_rpt_time
	null_sleep_stage_classify_acc_tmp = zeros(null_rept, 1);
    feature_data = feature_data_balance_cell{balance_i};
    for rpti = 1:null_rept
        disp(['null-model svm predict. randi:' num2str(rpti)]); tic;
        feature_temp = feature_data;
        y = feature_data(:, network_n+1);
        feature_temp(:, network_n+1) = y(randperm(length(y)));
        
        [externalAccuracies, externalConfusionMatrix, beta_Matrix, bias_Matrix, ~, ~, ~] = svm_predict_stage(feature_temp, 1, kfold, '');
        
        null_sleep_stage_classify_acc_tmp(rpti) = externalAccuracies;
        toc;
    end
    null_sleep_stage_classify_acc = [null_sleep_stage_classify_acc; null_sleep_stage_classify_acc_tmp];
end



save(['/nd_disk3/guoyuan/Jinlong/sleep_ret/HFR_ret/SVM_multiclass_control_' mask_name_list_w{mask_idx} '_50_rpt_classify_network' num2str(network_n) '.mat' ], ...
 	'sleep_stage_classify_acc', 'sleep_stage_classify_confMatrix', ...
     'sleep_stage_classify_beta', 'sleep_stage_classify_bias', ...
 	'null_sleep_stage_classify_acc', ...
    'externalConfusionMatrix_list', 'ROC_x', 'ROC_y', 'AUC_list');


function [externalAccuracies, externalConfusionMatrix, beta_Matrix, bias_Matrix, allTrueLabels, allPredictedScores, nc_time] = svm_predict_stage(data, nRepeat, kfold_time, ROC_file_stem)
	% Assuming X is the feature matrix and y is the label vector
    n = size(data, 2);
    X = data(:, 1:n-1); 
    y = data(:, n); 
    y_class = unique(y);
    y_class_n = length(y_class);

    
    k = kfold_time;  

    % Grid search hyperparameter range (e.g. C and gamma parameters of SVM)
    BoxConstraints = [0.1, 1, 10];  
    KernelScales = [0.1, 1, 10];    

    % Store the results of external cross validation
    externalAccuracies = zeros(nRepeat, 1);
    externalConfusionMatrix = zeros(y_class_n, y_class_n);
    allTrueLabels = [];  
    allPredictedScores = [];
    beta_Matrix = zeros(n-1, y_class_n, y_class_n);
    bias_Matrix = zeros(y_class_n, y_class_n);
    nc_time = 0;

    for repeat = 1:nRepeat
        % External cross validation partition
        cv = cvpartition(length(y), 'KFold', k);
        
        % Accuracy of internal cross validation storage (for hyperparameter tuning)
        foldAccuracies = zeros(cv.NumTestSets, 1);
        
        for fold = 1:cv.NumTestSets
            % Retrieve the indexes of the training and testing sets
            contain_all_class_flg = 100;
            while contain_all_class_flg
                trainIdx = cv.training(fold);
                testIdx = cv.test(fold);
                
                % Obtain training and testing set data
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
            
            % One-against-one Strategy: Train all possible binary classification SVMs
            classes = y_class;  
            numClasses = length(classes);
            models = cell(numClasses, numClasses);  
            
            % Train SVM models for all category pairs
            for i = 1:numClasses-1
                for j = i+1:numClasses
                    % Obtain training data for the current category pair
                    idx = (yTrain == classes(i)) | (yTrain == classes(j));
                    XPair = XTrain(idx, :);
                    yPair = yTrain(idx);
                    
                    % Convert the label to a binary label (class i is 1, class j is -1)
                    yPair(yPair == classes(i)) = 1;
                    yPair(yPair == classes(j)) = -1;
                    
                    % Training SVM model
                    bestSVMModel = svm_model_optimal(XPair, yPair, BoxConstraints, KernelScales);
                    models{i,j} = bestSVMModel;
                    beta_Matrix(:, i, j) = beta_Matrix(:, i, j) + bestSVMModel.Beta(:);
                    bias_Matrix(i, j) = bias_Matrix(i, j) + bestSVMModel.Bias;
                end
            end
            
            % Make predictions on the test set
            yPred = zeros(length(yTest), 1);
            decisionValues = zeros(length(yTest), numClasses);  % Store decision values for each category for each sample
            for i = 1:length(yTest)
                % Calculate which category each sample belongs to
                votes = zeros(numClasses, 1);
                for iClass = 1:numClasses-1
                    for jClass = iClass+1:numClasses
                        % Obtain the SVM model for the current category pair
                        model = models{iClass,jClass};
                        
                        % Predict the current test sample
                        testSample = XTest(i, :);
                        prediction = predict(model, testSample);
                        % Update decision values (decision values for each category)
                        if prediction == 1
                            decisionValues(i, iClass) = decisionValues(i, iClass) + 1;
                        else
                            decisionValues(i, jClass) = decisionValues(i, jClass) + 1;
                        end
                        
                        % Increase voting based on SVM prediction results
                        if prediction == 1
                            votes(iClass) = votes(iClass) + 1;
                        else
                            votes(jClass) = votes(jClass) + 1;
                        end
                    end
                end
                
                % Select the category with the most votes
                [~, maxVoteClass] = max(votes);
                yPred(i) = classes(maxVoteClass);
            end

            % Add real labels and predicted decision values to the storage list
            allTrueLabels = [allTrueLabels; yTest];
            allPredictedScores = [allPredictedScores; decisionValues];
            
            % Calculate the accuracy of this cross validation
            foldAccuracies(fold) = sum(yPred == yTest) / length(yTest);
            confMat = confusionmat(yTest, yPred);
            
            externalConfusionMatrix = externalConfusionMatrix + confMat;
            
            
            nc_time = nc_time + 1;
        end
        
        % Calculate the accuracy of external cross validation
        externalAccuracies(repeat) = mean(foldAccuracies);
    end

    % Calculate the final average accuracy and standard deviation
    avgAccuracy = mean(externalAccuracies);
    stdAccuracy = std(externalAccuracies);

    beta_Matrix = beta_Matrix ./ nc_time;
    bias_Matrix = bias_Matrix ./ nc_time;

    fprintf('average accuracy: %.4f\n', avgAccuracy);
    fprintf('accuracy standard deviation: %.4f\n', stdAccuracy);

    allPredictedScores = allPredictedScores ./ 6;
    if ~isempty(ROC_file_stem)
        draw_ROC(allTrueLabels, allPredictedScores, {'Awake', 'N1', 'N2', 'N3', 'REM'}, ROC_file_stem);
    end

end

function bestSVMModel = svm_model_optimal(XTrain, yTrain, BoxConstraints, KernelScales)
    % Perform grid search hyperparameter tuning
    bestAccuracy = -Inf;
    bestBoxConstraint = NaN;
    bestKernelScale = NaN;
    
    % Traverse each hyperparameter combination
    for boxC = BoxConstraints
        for kernelS = KernelScales
            % Train SVM model using current BoxConstraint and KernelScale
            SVMModel = fitcsvm(XTrain, yTrain, 'KernelFunction', 'linear', ...
                'BoxConstraint', boxC, 'KernelScale', kernelS, 'CrossVal', 'on');
            
            % Obtain the cross validation accuracy of the current model
            cvAccuracy = 1 - kfoldLoss(SVMModel);  % Calculate accuracy
            
            % Record the optimal hyperparameters
            if cvAccuracy > bestAccuracy
                bestAccuracy = cvAccuracy;
                bestBoxConstraint = boxC;
                bestKernelScale = kernelS;
            end
        end
    end
    
    % Train the final model using the optimal hyperparameters
    bestSVMModel = fitcsvm(XTrain, yTrain, 'KernelFunction', 'linear', ...
        'BoxConstraint', bestBoxConstraint, 'KernelScale', bestKernelScale);
end




function draw_ROC(allTrueLabels, allPredictedScores, stage_Names, ROC_file_stem)
    % Draw ROC curve
    fig = figure;
    hold on;
    y_class_n = length(unique(allTrueLabels));
    for classIdx = 1:y_class_n
        % Calculate the ROC curve for each category
        [X, Y, ~, AUC] = perfcurve(allTrueLabels, allPredictedScores(:, classIdx), classIdx, 'XCrit', 'fpr', 'YCrit', 'tpr');
        plot(X, Y, 'LineWidth', 2, 'DisplayName', [stage_Names{classIdx} ' (AUC = ' num2str(AUC) ')']);
    end
    xlabel('False Positive Rate');
    ylabel('True Positive Rate');
    title('ROC Curves for Multiclass');
    legend('Location', 'southeast');  
    grid on;
    hold off;
%     saveas(gcf, [ROC_file_stem '.png']);
    export_fig([ROC_file_stem '.png'], '-m6', '-q100');
    close(fig);
end
