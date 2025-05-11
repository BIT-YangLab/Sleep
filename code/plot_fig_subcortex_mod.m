function plot_fig_subcortex_mod(file_path, varargin)
% Visualizing subcortex
% Author: Xinyu Wu and Guoyuan Yang, BIT.
% 
% Inputs:
%     - file_path
%       Path of files to be plot, must be nifti format files with MNI152 2mm standard voxel space(109*91*91).
% 
%     - output_plot_path
%       Path to save figures, we recommend to save a png file.
%       It depends on the project export_fig.
%       about export_fig, please see: https://github.com/altmany/export_fig
% 
%     - color_map
%       A n*3 matrix contains RGB color.
%       This project has already contains project slanCM which contains numerous colormaps.
%       about slanCM, please see: https://github.com/slandarer/slanColor
% 
%     - max_scale
%       upper threshold to control the color to be plot.
% 
%     - min_scale
%       lower threshold to control the color to be plot.
%       
%     - structure_label
%       label number to be showed, using labels defined in the FreeSurfer and contained in 'Atlas_ROIs.2.nii.gz'.
%       structure label could be found in:
%           https://surfer.nmr.mgh.harvard.edu/fswiki/FsTutorial/AnatomicalROI/FreeSurferColorLUT.
%       useful infomations:
%           thalamus_label = [10, 49];
%           hippocampus_label = [17, 53];
%           cerebellum_label = [8, 47];
%           striatum_label = [11, 12, 13, 26, 50, 51, 52, 58];
%           all_subcortex_label = [8, 10, 11, 12, 13, 16, 17, 18, 26, 28, 47, 49, 50, 51, 52, 53, 54, 58, 60];

% TODO: 
% 1 - cmap and cscale.
% 2 - more flexible selection of the roi to plot. Could choose to use a nifti version ROI file or HCP grayordinate brainstructure labels.
%     if NO specify, plot every points which is not equal to zero.
% 3 - interpolate?
%     somebody asked me if i could add some interpolation while plotting the dots. why not be one option?
% Parse input
[basepath, ~, ~] = fileparts(which('plot_fig_subcortex'));
addpath(fullfile(basepath, 'dependencies', 'export_fig'));

p = inputParser;
addParameter(p, 'cmap', gray(256));
addParameter(p, 'cscale', [-1, 1]);
addParameter(p, 'nifti_roi_file', fullfile(basepath, 'rois', 'Atlas_ROIs.2.nii.gz'));
addParameter(p, 'nifti_roi_label', nan);
addParameter(p, 'export_fig_path', nan);
addParameter(p, 'interp', 'none');

parse(p, varargin{:});

min_scale = p.Results.cscale(1);
max_scale = p.Results.cscale(2);

file = MRIread(file_path);

rois = MRIread(p.Results.nifti_roi_file);

if isnan(p.Results.nifti_roi_label)
    rois.vol(rois.vol~=0) = 1;
else
    tmp_pos = ismember(rois.vol, p.Results.nifti_roi_label);
    rois.vol(tmp_pos) = 1;
    rois.vol(~tmp_pos) = 0;
end

if p.Results.interp ~= 'none'
    % interpolation, testing version
    [len_y, len_x, len_z] = size(rois.vol);
    [corase_x, corase_y, corase_z] = meshgrid(1:len_x, 1:len_y, 1:len_z);
    [fine_x, fine_y, fine_z] = meshgrid(1:0.5:len_x, 1:0.5:len_y, 1:0.5:len_z);
    interp_roi = interp3(corase_x, corase_y, corase_z, rois.vol, fine_x, fine_y, fine_z, 'nearest', 0);
    interp_data = interp3(corase_x, corase_y, corase_z, file.vol, fine_x, fine_y, fine_z, 'spline', 0);
    rois.vol = interp_roi;
    file.vol = interp_data;
end

roi_pos = find(rois.vol);

posx = zeros(length(roi_pos), 1);
posy = zeros(length(roi_pos), 1);
posz = zeros(length(roi_pos), 1);
color = zeros(length(roi_pos), 3);

scatter_size = 50;

for i = 1: length(roi_pos)
    [tmpx, tmpy, tmpz] = ind2sub([109, 91, 91], roi_pos(i));
    posx(i) = tmpx; posy(i) = tmpy; posz(i) = tmpz;

    if file.vol(roi_pos(i)) <= max_scale && file.vol(roi_pos(i)) > min_scale
        val = file.vol(roi_pos(i));
    elseif file.vol(roi_pos(i)) <= min_scale
        val = min_scale + (max_scale-min_scale)/size(color, 1);
    elseif file.vol(roi_pos(i)) > max_scale
        val = max_scale;
    end
        
    color(i, :) = p.Results.cmap(val, :);
end

figure;
scatter3(posx, posy, posz, scatter_size, color, 'filled');

% Adjusting figure to get better visualization
axis equal; axis vis3d;
% Adjusting figure size.
set(gcf, 'Position', [0, 0, 800, 800]);
% Hidding the axis.
set(gca, 'xticklabel', []); set(gca, 'yticklabel', []); set(gca, 'zticklabel', []);
set(gca, 'xtick', [], 'ytick', [], 'ztick', [], 'xcolor', 'w', 'ycolor', 'w', 'zcolor', 'w', 'Visible', 0);
% Setting background transparent.
set(gca, 'color', 'none'); set(gcf, 'color', 'none');

% Moving the camera position and angel to get better visualization.
% Recommendations:
%   for hippocampus and thalamus:
%       camorbit(95, -8, 'data', [0, 0, 1]);
%   for striatum:
%       camorbit(102, -2, 'data', [0, 0, 1]);
camorbit(102, -2, 'data', [0, 0, 1]);

% Detailed usage of export_fig, please see: https://github.com/altmany/export_fig
% '-m4' and '-q100' are used for sufficient figure quality.
% for better quality, try '-m10' and '-q100', but it would be slow to generate.
if ~isnan(p.Results.export_fig_path)
    export_fig(p.Results.export_fig_path, '-m4', '-q100');
end

end