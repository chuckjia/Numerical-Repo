% Graph_Soln.m - Graph numerical solutions and/or errors
%
% This script graphs the numerical solutions and errors from the humidity
% calculation.
%
% Author: Chuck Jia
% Created on: Oct 20, 2017

clear; clc

% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== %
% Plot Settings
% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== %

plotBoundary = false;
filename = "w_soln_32000";
titleText = "w at t = 10000.0s";


% contourLevels = [284, 290, 294, 296, 298, 300, 302, 308];  % For T
% contourLevels = [0.070, 0.883, 1.07, 1.25, 1.35, 1.45 1.55, 1.6, 1.7, 1.8] * 0.01;  % For q
% contourLevels = [5.6, 5.7, 5.8, 6, 6.5, 7, 8, 9, 10];
contourLevels = [-0.3, -0.025, -0.2, -0.1, 0, 0.1, 0.2, 0.025, 0.3];

% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== %
% Graph Selected Plots
% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== %

getPar_scp;  % Read and calculate parameters
getCellCenters_scp;  % Read cell centers
matShape = [numCellsX, numCellsP];
% Choose whether to plot the boundary or not
if (~plotBoundary)
    cellCentersX = rmBDVal_fcn(cellCentersX);
    cellCentersP = rmBDVal_fcn(cellCentersP);
end
fprintf("Plotting\n");

% Read numerical solutions/errors from file
folder = "MovieFrames/";
res = reshape(getVecFromFile_fcn(folder, filename), matShape);
if (~plotBoundary)
    res = rmBDVal_fcn(res);
end
% Plot the solution and error
near_mount_portion = 20 / 25;
viewAngle = [0, -90];
p_floor = floor(numCellsP * near_mount_portion); p_ceiling = numCellsP - 2 - 3;
x_coords = cellCentersX(p_floor:p_ceiling,:);
p_coords = cellCentersP(p_floor:p_ceiling,:);
res = res(p_floor:p_ceiling,:);
figure
surf(x_coords, p_coords, res);
y_axis_range = [650, pB];
xlim([x0, xf]); ylim(y_axis_range);
view(viewAngle)
figure
contourf(x_coords, p_coords, res, contourLevels);
xlim([x0, xf]); ylim(y_axis_range);
% colormap(hot)
view(viewAngle)
title(titleText)
% Add titles and labels

% title({graphTitleList(ff), titleLine2, ""});
xlabel('x coordinates'); ylabel('p coordinates');


clear titleLine2 matShape fileList graphTitleList graphList
