% Graph_Movie.m - Make movies for numerical solutions and/or errors
% 
% This script graphs the numerical solutions and errors from the humidity
% calculation and collects frames to make a movie.
%
% Author: Chuck Jia
% Created on: Oct 20, 2017

clear; clc

% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== %
% Plot Settings
% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== %

solnNo = 1;  % 1: T, 2: q, 3: u, 4:w
graphSelection = 2;  % 1: Numerical soln, 2: Error, 3: Exact soln
removeBoundaryVal = true;
plotFreq = 2;  % Choose the frequency of the frames

% Choose to manually set viewing angles
setViewingAngle = true; 
viewingAngle = [1 1 2];
%viewingAngle = 2;  % Flat; viewing from above

zAxisLimits = [-0.2 0.2];
figureSize = [10, 10, 900, 700];  % x-pos, y-pos, width, height

% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== %
% Plot The Soln/Error and Compile A Movie
% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== %

getPar_sct;  % Read and calculate parameters
getCellCenters_sct;  % Read cell centers
% File list and title settings
fileList = ["T_soln", "q_soln", "u_soln", "w_soln", ...
    "T_err", "q_err", "u_err", "w_err",...
    "T_exact", "q_exact", "u_exact", "w_exact" ];
graphTitleList = [
    "Numerical Solution: T", "Numerical Solution: q", ...
    "Numerical Solution: u", "Numerical Solution: w", ...
    "Error: T", "Error: q", "Error: u", "Error: w", ...
    "Exact Solution: T", "Exact Solution: q", ... 
    "Exact Solution: u", "Exact Solution: w"];
% Create list of graphs to be plotted
fileNo = solnNo + (graphSelection - 1) * 4;
filenamePrefix = fileList(fileNo);
graphTitle = graphTitleList(fileNo);
folder = "MovieFrames/"; 
matShape = [numCellsX, numCellsP];
% Choose whether to plot boundary values or not
if (removeBoundaryVal)
    cellCentersX = rmBDVal_fcn(cellCentersX);
    cellCentersP = rmBDVal_fcn(cellCentersP);
end
% Fix figure size
figure('pos', figureSize)
axisLimits = [x0 xf pA pBx0 zAxisLimits];
%numTimeSteps = 2;
frameNo = 0;

for tt = 1:numTimeSteps
    if mod(tt, plotFreq)
        continue
    end
    
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
    axis(axisLimits);
    if (setViewingAngle)
        view(viewingAngle);
    end
    caxis(axisLimits(5:6))
    % Add titles and labels
    title({graphTitle, "t = " + num2str(t) + "s"});
    xlabel('x coordinates'); ylabel('p coordinates');
    frameNo = frameNo + 1;
    F(frameNo) = getframe;
end

% Write frames to a movie
v = VideoWriter('Movies/T_soln.avi');
open(v)
writeVideo(v, F)
close(v)
