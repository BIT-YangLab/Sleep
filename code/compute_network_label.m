function network_all = compute_network_label(network_s, network_n)
% construct network label
% Author: Jinlong Li and Guoyuan Yang, BIT.
% 
% Inputs:
%     - network_s
%       The attribution mask of each network. vertex_n * network_n
% 
%     - network_n
%       number of network list
% Outputs:
%     - network_all
%       network label result 
    network_all = zeros(size(network_s, 1), 1);
    for ki = 1:network_n
        idx = find(network_s(:, ki) == 1);
        network_all(idx) = ki; 
    end
end