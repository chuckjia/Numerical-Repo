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

makeMovie = true;
movieName = 'T_err_normalAngle_100x100';

solnNo = 1;  % 1: T, 2: q, 3: u, 4:w
graphSelection = 2;  % 1: Numerical soln, 2: Error, 3: Exact soln
removeBoundaryVal = true;
plotFreq = 1;  % Choose the frequency of the frames

% Choose to manually set viewing angles
setViewingAngle = true;
viewingAngle = [1 1 2];  % Normal angle
%viewingAngle = [0 1 0];  % x side
%viewingAngle = [1 0 0];  % p side
%viewingAngle = 2;  % Flat; viewing from direct above

zAxisLimits = [-5e-5, 5e-5];
figureSize = [10, 10, 900, 700];  % x-pos, y-pos, width, height

numTimeStepLimit = 1e10;  % Limiting the number of time steps in plotting

% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== %
% Plot The Soln/Error and Compile A Movie
% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== %

getPar_sct;  % Read and calculate parameters
getCellCenters_sct;  % Read cell centers
% File list and title settings
graphLists_sct;
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
figure('pos', figureSize);
axisLimits = [x0 xf pA pBx0 zAxisLimits];
numTimeStepsToGraph = min(numTimeSteps, numTimeStepLimit);
frameNo = 0;
titleLine2 = "Mesh Size: " + int2str(Nx) + "x" +  int2str(Np) + ", " + ...
    "Dt = " + num2str(Dt) + "s";
for tt = 1:numTimeStepsToGraph
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
    title({graphTitle, titleLine2, "t = " + num2str(t) + "s", ""});
    xlabel('x coordinates'); ylabel('p coordinates');
    frameNo = frameNo + 1;
    if (makeMovie)
        F(frameNo) = getframe(gcf);
        if (tt == 1) % Allocate memory space at the beginning
            numFrames = floor(numTimeStepsToGraph / plotFreq);
            singleFrame = F(frameNo);
            F = repmat(singleFrame, 1, numFrames);
        end
    else
        drawnow update;
    end
end

% Write frames to a movie
if (makeMovie)
    v = VideoWriter(strcat('Movies/', movieName));
    open(v);
    writeVideo(v, F);
    close(v);
end
 clear axisLimits F figureSize fileList filename filenamePrefix fileNo ...
    folder frameNo graphSelection graphTitle graphTitleList makeMovie ...
    matShape movieName numTimeStepLimit numTimeStepsToGraph plotFreq ...
    removeBoundaryVal res singleFrame setViewingAngle numFrames solnNo ... 
    t titleLine2 tt v viewingAngle zAxisLimits