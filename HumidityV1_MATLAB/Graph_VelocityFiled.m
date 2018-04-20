% Graph_VelocityField.m - Grap velocity fields
%
% This script graphs the 
%
% Author: Chuck Jia
% Created on: Nov 15, 2017

clear; clc

% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== %
% Plot Settings
% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== %

plot_time = 6000
graphSelection = 1;  % 1: Numerical soln, 2: Exact
removeBoundaryVal = false;

%viewingAngle = [1 1 2];  % Normal angle
viewingAngle = 2;  % Flat; viewing from direct above

figureSize = [10, 10, 900, 700];  % x-pos, y-pos, width, height

% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== %
% Plot The Soln/Error and Compile Movie
% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== %

folder_getParScp = "/Users/chuckjia/Documents/Workspace/DataStorage/Humidity/Sim10/Norm/";
getPar_scp;  % Read and calculate parameters
getCellCenters_scp;  % Read cell centers

% Create list of graphs to be plotted
if (graphSelection == 1)
    filenamePrefix = ["u_soln", "w_soln"];
    graphTitle = 'Numerical Velocity Field';
else
    filenamePrefix = ["u_exact", "w_exact"];
    graphTitle = 'Exact Velocity Field';
end
folder = "/Users/chuckjia/Documents/Workspace/DataStorage/Humidity/Sim09/Snapshots/";
matShape = [numCellsX, numCellsP];
quiverScaleOnX = 0.05;  % This is introduced to adjust arrow head angle
cellCentersX = cellCentersX .* quiverScaleOnX;
% Choose whether to plot boundary values or not
if (removeBoundaryVal)
    cellCentersX = rmBDVal_fcn(cellCentersX);
    cellCentersP = rmBDVal_fcn(cellCentersP);
end
% Fix figure size
figure('pos', figureSize);
axisLimits = [x0, xf * quiverScaleOnX, pA, pBx0];
titleLine2 = "Mesh Size: " + int2str(Nx) + "x" +  int2str(Np) + ", " + ...
    "Time = " + num2str(plot_time) + "s";


filename = "u_soln_" + num2str(plot_time * 2)
% Read numerical solutions/errors from file
u_res = reshape(getVecFromFile_fcn(folder, filename), matShape);
filename = "w_soln_" + num2str(plot_time * 2)
w_res = reshape(getVecFromFile_fcn(folder, filename), matShape);
if (removeBoundaryVal)
    u_res = rmBDVal_fcn(u_res);
    w_res = rmBDVal_fcn(w_res);
end
% Plot the solution and error
xRange = 1:10; pRange = 180:200; 
h = quiver(cellCentersX(pRange, xRange), cellCentersP(pRange, xRange), u_res(pRange, xRange), w_res(pRange, xRange));
% axis(axisLimits);
view(viewingAngle);
% Add titles and labels
title({graphTitle, titleLine2});
xlabel('x coordinates'); ylabel('p coordinates');

% Print graph to pdf with the best fit page size
print('velocityField_numer','-dpdf','-bestfit');  


clear axisLimits cellCentersP cellCentersX figureSize filename ...
    filenamePrefix folder graphSelection graphTitle h matShape ... 
    quiverScaleOnX removeBoundaryVal titleLine2 u_res w_res
