function [dice_list, label_list] = compute_dice_coef(label1, label2)
% compute dice coefficients
% Author: Jinlong Li and Guoyuan Yang, BIT.
% 
% Inputs:
%     - label1, label2
%       Two atlas for calculating Dice coefficient (same dimension)
% 
% Outputs:
%     - dice_list
%       list of dice coefficient
%     - label_list
%       list of network or parcellations
    label_list = union(unique(label1), unique(label2));
    label_list(label_list < 1) = [];
    label_list(isnan(label_list)) = [];

    dice_list = zeros(length(label_list), 1);
    for pi = 1:length(label_list)
        index1 = label1 == label_list(pi);
        index2 = label2 == label_list(pi);
        index = index1 + index2;
        dice_list(pi) = nnz(index == 2) / nnz(index);
    end
end