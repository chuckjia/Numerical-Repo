% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
% Graph Numerical Solution or Errors
% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== %
clear; clc

% Read and calculate parameters
getPar_sct;
% Read cell centers
getCellCenters_sct;
% File list and title settings
fileList = ["T_soln", "q_soln", "u_soln", "w_soln", ...
    "T_err", "q_err", "u_err", "w_err",...
    "T_exact", "q_exact", "u_exact", "w_exact" ];
graphTitles = [
    "Numerical Solution: T", "Numerical Solution: q", ...
    "Numerical Solution: u", "Numerical Solution: w", ...
    "Error: T", "Error: q", "Error: u", "Error: w", ...
    "Exact Solution: T", "Exact Solution: q", "Exact Solution: u", "Exact Solution: w"];
numGraphs = size(fileList);

% Choose graphs to plot
graphNumericalSoln = [1];
graphErr = [1];
graphExactSoln = [];
removeBoundaryVal = true;

% List of graphs to be plotted
graphList = [graphNumericalSoln, graphErr + 4, graphExactSoln + 8];
folder = "Results/"; matShape = [numCellsX, numCellsP];
if (removeBoundaryVal)
    cellCentersX = rmBDVal_fcn(cellCentersX);
    cellCentersP = rmBDVal_fcn(cellCentersP);
end
for ff = graphList
    % Read numerical solutions/errors from file
    res = reshape(getVecFromFile_fcn(folder, fileList(ff)), matShape);
    if (removeBoundaryVal)
        res = rmBDVal_fcn(res);
    end
    % Plot the solution and error
    figure
    surf(cellCentersX, cellCentersP, res);
    % Add titles and labels
    title({graphTitles(ff)});
    xlabel('x coordinates'); ylabel('p coordinates');
end
