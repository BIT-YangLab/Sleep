% Author: Jinlong Li and Guoyuan Yang, BIT.
% classification for sleep stage with network size ratio by ANN
% 
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

load([ret_dir '/IndiPar_net_multi_7/Network_ind_subcortex_cleanup8min.mat']);
addpath('/nd_disk3/guoyuan/Xinyu/plot_fig_subcortex/dependencies/export_fig/');

network_n = 7;
stageNames = {'Awake', 'N1', 'N2', 'N3', 'REM'};
rept_svm = 40;
kfold = 5;
nTrees = 100;
null_rept = 50;
balance_rpt_time = 5;
hiddenLayerSize = 10; 

% create match table
% network_label_all = network_label_subcortex_lisa;
network_label_all = network_label_subcortex_lisa;
feature_data = [];

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
load([ret_dir '/DICE_contral_7.mat']);


for stagei = 1:4
    for subi = 1:size(network_label_all, 2)
        if isempty(network_label_all{stagei, subi})
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

for balance_i = 1:balance_rpt_time
    feature_data = feature_data_balance_cell{balance_i};
    [externalAccuracies, externalConfusionMatrix, allTrueLabels, allPredictedScores] = ANN_predict_stage(feature_data, rept_svm, kfold, hiddenLayerSize, ...
    '');
    sleep_stage_classify_acc = [sleep_stage_classify_acc; externalAccuracies];
    sleep_stage_classify_confMatrix{balance_i} = externalConfusionMatrix;

    allTrueLabels_list = [allTrueLabels_list; allTrueLabels];
    allPredictedScores_list = [allPredictedScores_list; allPredictedScores];

    %draw_confMatrix(externalConfusionMatrix, stageNames, ...
    %    ['/nd_disk3/guoyuan/Jinlong/sleep_ret/HFR_ret/pic/SVM_control_cbm_Mean_ConfMatrix_heatmap_Multiclass_network' num2str(network_n) '_balance' num2str(balance_i) '.png']);
		
end

% ROC curve 
['/nd_disk3/guoyuan/Jinlong/sleep_ret/HFR_ret/pic/SVM_control_cbm_' num2str(network_n) '_ROC_balance' num2str(balance_i)]
y_class_n = length(unique(allTrueLabels_list));
ROC_x = [];
ROC_y = [];
AUC_list = [];
for classIdx = 1:y_class_n
    % compute ROC curve
    [X, Y, ~, AUC] = perfcurve(allTrueLabels_list, allPredictedScores_list(:, classIdx), classIdx, 'XCrit', 'fpr', 'YCrit', 'tpr');
    
    % downsample
    indices = round(linspace(1, length(X), 100));
    
    ROC_x_reduced = X(indices);
    ROC_y_reduced = Y(indices);
    

    ROC_x = [ROC_x; ROC_x_reduced];
    ROC_y = [ROC_y; ROC_y_reduced];
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
        
        [externalAccuracies, externalConfusionMatrix, ~, ~] = ANN_predict_stage(feature_temp, 1, kfold, hiddenLayerSize, '');
        
        null_sleep_stage_classify_acc_tmp(rpti) = externalAccuracies;
        toc;
    end
    null_sleep_stage_classify_acc = [null_sleep_stage_classify_acc; null_sleep_stage_classify_acc_tmp];
end


% leave one network out prediction
leave1_sleep_stage_classify_acc = [];
leave1_sleep_stage_classify_acc_org = [];
leave2_sleep_stage_classify_acc = [];
leave1_sleep_stage_classify_confMatrix = cell(network_n, balance_rpt_time);

save(['/nd_disk3/guoyuan/Jinlong/sleep_ret/HFR_ret/ANN_multiclass_' mask_name_list_w{mask_idx} '_50_rpt_classify_network' num2str(network_n) '.mat' ], ...
 	'sleep_stage_classify_acc', 'sleep_stage_classify_confMatrix', ...
     'sleep_stage_classify_beta', 'sleep_stage_classify_bias', ...
 	'null_sleep_stage_classify_acc', ...
 	'leave1_sleep_stage_classify_acc',  'leave1_sleep_stage_classify_confMatrix', ...
    'ROC_x', 'ROC_y', 'AUC_list');

function [externalAccuracies, externalConfusionMatrix, allTrueLabels, allPredictedScores] = ANN_predict_stage(data, nRepeat, kfold_time, hiddenLayerSize, ROC_file_stem)
	% X: feature, Y:label
    n = size(data, 2);
    X = data(:, 1:n-1); % feature
    y = data(:, n); % label
    Y = full(ind2vec(y'))';  % one-hot code
    y_class = unique(y);
    y_class_n = length(y_class);

    % k-fold
    k = kfold_time; 

    % save results
    externalAccuracies = zeros(nRepeat, 1);
    externalConfusionMatrix = zeros(y_class_n, y_class_n);
    allTrueLabels = [];  % true label
    allPredictedScores = [];  % predicted socres
    nc_time = 0;

    for repeat = 1:nRepeat
        % external cross validation
        cv = cvpartition(length(y), 'KFold', k);
        
        % Accuracy of internal cross validation storage (for hyperparameter tuning)
        foldAccuracies = zeros(cv.NumTestSets, 1);

        for i = 1:cv.NumTestSets
            % Divide the training set and testing set
            trainIdx = cv.training(i);
            testIdx = cv.test(i);
            
            % Training a neural network using a training set
            net = feedforwardnet(hiddenLayerSize);
            net.trainFcn = 'trainlm';
            net.layers{end}.transferFcn = 'softmax';
            net = configure(net, X(trainIdx, :)', Y(trainIdx, :)');
            [net, ~] = train(net, X(trainIdx, :)', Y(trainIdx, :)');
            
            % Predictive testing set
            output = net(X(testIdx, :)');
            [~, predictedLabels] = max(output, [], 1);
            predictedLabels = predictedLabels';
            
            % Calculate accuracy
            foldAccuracies(i) = sum(predictedLabels == vec2ind(Y(testIdx, :)')') / length(Y(testIdx));

            confMat = confusionmat(vec2ind(Y(testIdx, :)'), predictedLabels);
            
            externalConfusionMatrix = externalConfusionMatrix + confMat;

            % Store all true labels and predicted probabilities (for plotting ROC curves)
            allTrueLabels = [allTrueLabels; vec2ind(Y(testIdx, :)')'];  % Real label
            allPredictedScores = [allPredictedScores; output'];  % Prediction probability for each category

            nc_time = nc_time + 1;
        end
        
        % Calculate the accuracy of external cross validation
        externalAccuracies(repeat) = mean(foldAccuracies);
    end

    % Calculate the final average accuracy and standard deviation
    avgAccuracy = mean(externalAccuracies);
    stdAccuracy = std(externalAccuracies);

    % confusion matrix 0 - 1
    externalConfusionMatrix = externalConfusionMatrix ./ nc_time;
    externalConfusionMatrix = externalConfusionMatrix ./ sum(externalConfusionMatrix, 2);


    fprintf('average accuracy: %.4f\n', avgAccuracy);
    fprintf('accuracy standard deviation: %.4f\n', stdAccuracy);

    if ~isempty(ROC_file_stem)
        draw_ROC(allTrueLabels, allPredictedScores, {'Awake', 'N1', 'N2', 'N3', 'REM'}, ROC_file_stem);
    end
    


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