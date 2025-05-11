% Author: Jinlong Li and Guoyuan Yang, BIT.
% binary classification for sleep stage with network size ratio by SVM
% 
clc;clear;

work_dir = '/nd_disk3/guoyuan/sleep/a_bash/Gordon_method/HFR_ai_mod/';
ret_dir = '/nd_disk3/guoyuan/Jinlong/sleep_ret/HFR_ret/';

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
stageNames = {'Awake', 'N1', 'N2', 'N3', 'REM'};
rept_svm = 40;
kfold = 5;
balance_rpt_time = 5;


% roi mask
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
mask_idx = 4;

% create match table
network_label_all = network_label_subcortex_lisa;
feature_data = [];

for stagei = 1:4
	for subi = 1:size(network_label_all, 2)
		if isempty(network_label_all{stagei, subi})
			continue;
		end
		for sessi = 1:length(network_label_all{stagei, subi})
			label_tmp = network_label_all{stagei, subi}{sessi};
			label_tmp(mask_list_w{mask_idx} == 0) = 0;
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
sleep_stage_classify_acc = cell(5, 5);
sleep_stage_classify_std = zeros(5, 5);
sleep_stage_classify_confMatrix = cell(5, 5);
sleep_stage_classify_beta = zeros(5, 5, network_n);
sleep_stage_classify_bias = zeros(5, 5);
for stagei = 1:3
	for stagej = stagei+1:4
        acc_list_tmp = [];
        for balance_i = 1:balance_rpt_time
            feature_data = feature_data_balance_cell{balance_i};
            disp(['svm predict - ' num2str(balance_i) ' : ' stageNames{stagei} ' vs ' stageNames{stagej}]); tic;
		    feature_temp = feature_data(:, 1:network_n);
		    label_bi = feature_data(:, network_n+1) == stagei;
    
		    % filter sleep stage
		    label_bj = feature_data(:, network_n+1) == stagej;
		    label_filter = (label_bi + label_bj) == 0;
		    feature_temp = [feature_temp label_bi];
		    feature_temp(label_filter, :) = [];
    
		    [externalAccuracies, externalConfusionMatrix, beta_Matrix, bias_Matrix] = svm_predict_stage(feature_temp, rept_svm, kfold);
		    acc_list_tmp = [acc_list_tmp; externalAccuracies];

		
		    toc;
        end
        sleep_stage_classify_acc{stagei, stagej} = acc_list_tmp;

	end
end

% null model prediction
null_rept = 50;
null_sleep_stage_classify_acc = zeros(5, 5, null_rept);
null_sleep_stage_classify_std = zeros(5, 5, null_rept);
null_sleep_stage_classify_confMatrix = cell(5, 5, null_rept);
for stagei = 1:3
	for stagej = stagei+1:4
        for balance_i = 1:balance_rpt_time
            tic;
            feature_data = feature_data_balance_cell{balance_i};
		    feature_temp = feature_data(:, 1:network_n);
		    label_bi = feature_data(:, network_n+1) == stagei;
    
		    % filter sleep stage
		    label_bj = feature_data(:, network_n+1) == stagej;
		    label_filter = (label_bi + label_bj) == 0;
		    feature_temp = [feature_temp label_bi];
		    feature_temp(label_filter, :) = [];
		    
		    % random sleep stage label
		    for rpti = 1:null_rept
			    disp(['null-model svm predict - ' num2str(balance_i) ': ' stageNames{stagei} ' vs ' stageNames{stagej} ' randi:' num2str(rpti)]); 
			    y = feature_temp(:, network_n+1);
			    feature_temp(:, network_n+1) = y(randperm(length(y)));
			    
			    [externalAccuracies, externalConfusionMatrix, beta_Matrix, bias_Matrix] = svm_predict_stage(feature_temp, 1, kfold);
			    
			    null_sleep_stage_classify_acc(stagei, stagej, rpti+(balance_i-1)*null_rept) = externalAccuracies;
			    
            end
            toc;
        end
	end
end



save(['/nd_disk3/guoyuan/Jinlong/sleep_ret/HFR_ret/SVM_40rpt_binary_' mask_name_list_w{mask_idx} '_classify_network' num2str(network_n) '.mat' ], ...
	'sleep_stage_classify_acc', 'sleep_stage_classify_confMatrix', ...
    'sleep_stage_classify_beta', 'sleep_stage_classify_bias', ...
	'null_sleep_stage_classify_acc', 'null_sleep_stage_classify_std', 'null_sleep_stage_classify_confMatrix');


function [externalAccuracies, externalConfusionMatrix, beta_Matrix, bias_Matrix] = svm_predict_stage(data, nRepeat, kfold_time)
	% Assuming X is the feature matrix and y is the label vector
	n = size(data, 2);
	X = data(:, 1:n-1); 
	y = data(:, n); 

	% Set external cross validation parameters
	k = kfold_time;  

	% Grid search hyperparameter range (e.g. C and gamma parameters of SVM)
	BoxConstraints = [0.1, 1, 10];  % Select candidate values for BoxConstraint
	KernelScales = [0.1, 1, 10];    % Select candidate values for KernelScale

	% Store the results of external cross validation
	externalAccuracies = zeros(nRepeat, 1);
	externalConfusionMatrix = zeros(2, 2);
	beta_Matrix = zeros(n-1, 1);
	bias_Matrix = 0;
	nc_time = 0;

	for repeat = 1:nRepeat
		% External cross validation partition
		cv = cvpartition(length(y), 'KFold', k);
		
		% Accuracy of internal cross validation storage (for hyperparameter tuning)
		foldAccuracies = zeros(cv.NumTestSets, 1);
		
		for fold = 1:cv.NumTestSets
			% Retrieve the indexes of the training and testing sets
			trainIdx = cv.training(fold);
			testIdx = cv.test(fold);
			
			% Obtain training and testing set data
			XTrain = X(trainIdx, :);
			yTrain = y(trainIdx);
			XTest = X(testIdx, :);
			yTest = y(testIdx);
			
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
					cvAccuracy = 1 - kfoldLoss(SVMModel);  
					
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
			
			% Evaluate the model on the test set
			yPred = predict(bestSVMModel, XTest);
			foldAccuracies(fold) = sum(yPred == yTest) / length(yTest);
			confMat = confusionmat(yTest, yPred);
			
			externalConfusionMatrix = externalConfusionMatrix + confMat;
			beta_Matrix = beta_Matrix + bestSVMModel.Beta;
			bias_Matrix = bias_Matrix + bestSVMModel.Bias;

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

	beta_Matrix = beta_Matrix ./ nc_time;
	bias_Matrix = bias_Matrix ./ nc_time;

	fprintf('average accuracy: %.4f\n', avgAccuracy);
	fprintf('accuracy standard deviation: %.4f\n', stdAccuracy);


end

