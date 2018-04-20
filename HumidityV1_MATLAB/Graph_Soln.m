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

graphNumerSoln = [1];
graphErr = [1];
graphExactSoln = [];
plotBoundary = false;

% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== %
% Graph Selected Plots
% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== %

folder = "/Users/chuckjia/Documents/Workspace/Git/Numerical-Repo/HumidityV1/"
getPar_scp;  % Read and calculate parameters
getCellCenters_scp;  % Read cell centers
% File list and title settings
graphLists_scp;
% List of graphs to be plotted
graphList = [graphNumerSoln, graphErr + 4, graphExactSoln + 8];
folder = "Results/"; 
matShape = [numCellsX, numCellsP];
% Choose whether to plot the boundary or not
if (~plotBoundary)
    cellCentersX = rmBDVal_fcn(cellCentersX);
    cellCentersP = rmBDVal_fcn(cellCentersP);
end
titleLine2 = "Mesh Size: " + int2str(Nx) + "x" +  int2str(Np) + ", " + ... 
    "Dt = " + num2str(Dt) + "s";
fprintf("Plotting\n");
for ff = graphList
    % Read numerical solutions/errors from file
    fprintf("  - " + graphTitleList(ff) + "\n");
    res = reshape(getVecFromFile_fcn(folder, fileList(ff)), matShape);
    if (~plotBoundary)
        res = rmBDVal_fcn(res);
    end
    % Plot the solution and error
    figure; surf(cellCentersX, cellCentersP, res);
    % Add titles and labels
    
    title({graphTitleList(ff), titleLine2, ""});
    xlabel('x coordinates'); ylabel('p coordinates');
end
%view(-180, 90)

clear titleLine2 matShape fileList graphTitleList graphList
