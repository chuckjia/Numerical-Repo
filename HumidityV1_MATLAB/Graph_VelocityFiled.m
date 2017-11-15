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

graphSelection = 1;  % 1: Numerical soln, 2: Exact
removeBoundaryVal = true;

%viewingAngle = [1 1 2];  % Normal angle
viewingAngle = 2;  % Flat; viewing from direct above

figureSize = [10, 10, 900, 700];  % x-pos, y-pos, width, height

% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== %
% Plot The Soln/Error and Compile Movie
% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== %

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
folder = "Results/";
matShape = [numCellsX, numCellsP];
quiverScaleOnX = 0.05;  % This is introduced to adjust arrow head angle
% Choose whether to plot boundary values or not
if (removeBoundaryVal)
    cellCentersX = rmBDVal_fcn(cellCentersX) .* quiverScaleOnX;
    cellCentersP = rmBDVal_fcn(cellCentersP);
end
% Fix figure size
figure('pos', figureSize);
axisLimits = [x0, xf * quiverScaleOnX, pA, pBx0];
titleLine2 = "Mesh Size: " + int2str(Nx) + "x" +  int2str(Np) + ", " + ...
    "Dt = " + num2str(Dt) + "s";

filename = filenamePrefix(1)
% Read numerical solutions/errors from file
u_res = reshape(getVecFromFile_fcn(folder, filename), matShape);
filename = filenamePrefix(2)
w_res = reshape(getVecFromFile_fcn(folder, filename), matShape);
if (removeBoundaryVal)
    u_res = rmBDVal_fcn(u_res);
    w_res = rmBDVal_fcn(w_res);
end
% Plot the solution and error
h = quiver(cellCentersX, cellCentersP, u_res, w_res);
axis(axisLimits);
view(viewingAngle);
% Add titles and labels
title({graphTitle, titleLine2});
xlabel('x coordinates'); ylabel('p coordinates');

% Print graph to pdf with the best fit page size
print('velocityField_numer','-dpdf','-bestfit');  


clear axisLimits cellCentersP cellCentersX figureSize filename ...
    filenamePrefix folder graphSelection graphTitle h matShape ... 
    quiverScaleOnX removeBoundaryVal titleLine2 u_res w_res
