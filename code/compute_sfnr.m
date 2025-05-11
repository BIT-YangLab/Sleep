function tmp_sfnr = compute_sfnr(tmp_indiv_data)
% Compute sfnr map for cifti file
% Author: Jinlong Li and Guoyuan Yang, BIT.
% 
% Inputs:
%     - tmp_indiv_data
%       timeseries vertex_n*timepoints
% Outputs:
%     - tmp_sfnr
%       sfnr map

    % tmp_indiv_data: vertex * timepoints
    mean_tmp_indiv_data = mean(tmp_indiv_data, 2);
    std_tmp_indiv_data = std(tmp_indiv_data, 0, 2);
    
    if mean_tmp_indiv_data == 0 
        if std_tmp_indiv_data == 0
            error('wrong data');
        end
    end
    
    tmp_sfnr = mean_tmp_indiv_data ./ std_tmp_indiv_data;

    tmp_sfnr(isnan(tmp_sfnr)) = 0;

end