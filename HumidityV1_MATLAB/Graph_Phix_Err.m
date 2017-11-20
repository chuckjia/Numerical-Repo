% Graph_Phix_Err.m - Graph numerical solutions and/or errors
%
% This script graphs the numerical solutions and errors from the humidity
% calculation.
%
% Author: Chuck Jia
% Created on: Nov 19, 2017

clear; clc

% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== %
% Graph Selected Plots
% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== %

getPar_scp;  % Read and calculate parameters
getCellCenters_scp;  % Read cell centers
folder = "Results/";
matShape = [Nx, Np];
cellCentersX = rmBDVal_fcn(cellCentersX);
cellCentersP = rmBDVal_fcn(cellCentersP);
titleLine2 = "Mesh Size: " + int2str(Nx) + "x" +  int2str(Np) + ", " + ...
    "Dt = " + num2str(Dt) + "s";
fprintf("Plotting\n");

% Read numerical solutions/errors from file
res = reshape(getVecFromFile_fcn(folder, "phix_err"), matShape);
% Plot the solution and error
figure; surf(cellCentersX, cellCentersP, res);
% Add titles and labels

title({"\phi_x Error", titleLine2, ""});
xlabel('x coordinates'); ylabel('p coordinates');


clear titleLine2 matShape fileList graphTitleList graphList
