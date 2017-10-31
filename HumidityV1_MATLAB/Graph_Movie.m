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
graphTitleList = [
    "Numerical Solution: T", "Numerical Solution: q", ...
    "Numerical Solution: u", "Numerical Solution: w", ...
    "Error: T", "Error: q", "Error: u", "Error: w", ...
    "Exact Solution: T", "Exact Solution: q", "Exact Solution: u", "Exact Solution: w"];
numGraphs = size(fileList);

% Choose graphs to plot
solnNo = 1;
graphSelection = 1;  % 1: Numerical solution, 2: Error, 3: Exact solution
removeBoundaryVal = true;

% List of graphs to be plotted
fileNo = solnNo + (graphSelection - 1) * 4;
filenamePrefix = fileList(fileNo);
graphTitle = graphTitleList(fileNo);
folder = "MovieFrames/"; matShape = [numCellsX, numCellsP];
if (removeBoundaryVal)
    cellCentersX = rmBDVal_fcn(cellCentersX);
    cellCentersP = rmBDVal_fcn(cellCentersP);
end
setViewingAngle = true;
viewingAngle = [1 1 2];
%viewingAngle = 2; % Flat

figure('pos', [10, 10, 900, 700])
%numTimeSteps = 2;
for tt = 1:numTimeSteps
    t = tt * Dt;
    fprintf("Currently graphing plot no. %d\n", tt);
    filename = filenamePrefix + "_" + int2str(tt);
    % Read numerical solutions/errors from file
    res = reshape(getVecFromFile_fcn(folder, filename), matShape);
    if (removeBoundaryVal)
        res = rmBDVal_fcn(res);
    end
    % Plot the solution and error
    surf(cellCentersX, cellCentersP, res);
    axis([0 50000 200 1000 -5 5]);
    if (setViewingAngle)
        view(viewingAngle);
    end
    caxis([-5, 5])
    % Add titles and labels
    title({graphTitle, "t = " + num2str(t) + "s"});
    xlabel('x coordinates'); ylabel('p coordinates');
    F(tt) = getframe;
end

v = VideoWriter('Movies/T_soln.avi');
open(v)
writeVideo(v, F)
close(v)
