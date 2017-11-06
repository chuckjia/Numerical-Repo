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
graphExactSoln = [1];
plotBoundary = false;

% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== %
% Graph Selected Plots
% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== %

getPar_sct;  % Read and calculate parameters
getCellCenters_sct;  % Read cell centers
% File list and title settings
graphLists_sct;
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
for ff = graphList
    % Read numerical solutions/errors from file
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

clear titleLine2 matShape fileList graphTitleList graphList