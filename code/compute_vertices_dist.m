function dist_W = compute_vertices_dist(roi_vertices_pos_i, roi_vertices_pos_j)
% compute distance between two vertex in fs_LR_32k
% Author: Jinlong Li and Guoyuan Yang, BIT.
% 
% Inputs:
%     - roi_vertices_pos_i/roi_vertices_pos_j:
%       index of vertex in fs_LR_32k
% 
% Outputs:
%     - dist_W
%       Euclidean distance

    dist_x = repmat(roi_vertices_pos_i(:, 1), 1, size(roi_vertices_pos_j, 1)) - repmat(roi_vertices_pos_j(:, 1), 1, size(roi_vertices_pos_i, 1))';
    dist_y = repmat(roi_vertices_pos_i(:, 2), 1, size(roi_vertices_pos_j, 1)) - repmat(roi_vertices_pos_j(:, 2), 1, size(roi_vertices_pos_i, 1))';
    dist_z = repmat(roi_vertices_pos_i(:, 3), 1, size(roi_vertices_pos_j, 1)) - repmat(roi_vertices_pos_j(:, 3), 1, size(roi_vertices_pos_i, 1))';

    dist_W = sqrt(dist_x.^2 + dist_y.^2 + dist_z.^2);

end