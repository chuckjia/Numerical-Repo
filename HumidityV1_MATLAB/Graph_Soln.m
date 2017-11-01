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
fileList = ["T_soln", "q_soln", "u_soln", "w_soln", ...
    "T_err", "q_err", "u_err", "w_err",...
    "T_exact", "q_exact", "u_exact", "w_exact" ];
graphTitles = [
    "Numerical Solution: T", "Numerical Solution: q", ...
    "Numerical Solution: u", "Numerical Solution: w", ...
    "Error: T", "Error: q", "Error: u", "Error: w", ...
    "Exact Solution: T", "Exact Solution: q", ... 
    "Exact Solution: u", "Exact Solution: w"];
% List of graphs to be plotted
graphList = [graphNumerSoln, graphErr + 4, graphExactSoln + 8];
folder = "Results/"; 
matShape = [numCellsX, numCellsP];
% Choose whether to plot the boundary or not
if (~plotBoundary)
    cellCentersX = rmBDVal_fcn(cellCentersX);
    cellCentersP = rmBDVal_fcn(cellCentersP);
end
for ff = graphList
    % Read numerical solutions/errors from file
    res = reshape(getVecFromFile_fcn(folder, fileList(ff)), matShape);
    if (~plotBoundary)
        res = rmBDVal_fcn(res);
    end
    % Plot the solution and error
    figure; surf(cellCentersX, cellCentersP, res);
    % Add titles and labels
    title({graphTitles(ff),""});
    xlabel('x coordinates'); ylabel('p coordinates');
end
