% Graph_VelocityField.m - Make movies for velocity fields
%
% This script graphs the velocity field from the humidity calculation and
% collects frames to make a movie.
%
% Author: Chuck Jia
% Created on: Oct 20, 2017

clear; clc

% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== %
% Plot Settings
% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== %

makeMovie = false;
movieName = 'velocityField_50x50';

graphSelection = 2;  % 1: Numerical soln, 2: Exact, 3: u, 4: w
removeBoundaryVal = true;
plotFreq = 1;  % Choose the frequency of the frames

% Choose to manually set viewing angles
setViewingAngle = true;
%viewingAngle = [1 1 2];  % Normal angle
viewingAngle = 2;  % Flat; viewing from direct above

% zAxisLimits = [-5e-5, 5e-5];
zAxisLimits = [-5, 5];
figureSize = [10, 10, 900, 700];  % x-pos, y-pos, width, height

numTimeStepLimit = 1;  % Limiting the number of time steps in plotting

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
folder = "MovieFrames/";
matShape = [numCellsX, numCellsP];
quiverScaleOnX = 0.05;  % This is introduced to adjust arrow head angle
% Choose whether to plot boundary values or not
if (removeBoundaryVal)
    cellCentersX = rmBDVal_fcn(cellCentersX) .* quiverScaleOnX;
    cellCentersP = rmBDVal_fcn(cellCentersP);
end
% Fix figure size
figure('pos', figureSize);
axisLimits = [x0, xf * quiverScaleOnX, pA, pBx0, zAxisLimits];
numTimeStepsToGraph = min(numTimeSteps, numTimeStepLimit);
frameNo = 0;
titleLine2 = "Mesh Size: " + int2str(Nx) + "x" +  int2str(Np) + ", " + ...
    "Dt = " + num2str(Dt) + "s";
cPlotCount = 0;
for tt = 1:numTimeStepsToGraph
    if mod(tt, cMovieFrameFreq)
        continue
    else
        cPlotCount = cPlotCount + 1;
    end
    if mod(cPlotCount, plotFreq)
        continue
    end
    
    t = tt * Dt;
    fprintf("Currently graphing plot no. %d\n", tt);
    filename = filenamePrefix(1) + "_" + int2str(tt);
    % Read numerical solutions/errors from file
    u_res = reshape(getVecFromFile_fcn(folder, filename), matShape);
    filename = filenamePrefix(2) + "_" + int2str(tt);
    w_res = reshape(getVecFromFile_fcn(folder, filename), matShape);
    if (removeBoundaryVal)
        u_res = rmBDVal_fcn(u_res);
        w_res = rmBDVal_fcn(w_res);
    end
    % Plot the solution and error
    h = quiver(cellCentersX, cellCentersP, u_res, w_res);
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

% print('velocityField','-dpdf','-bestfit');  // Print graph to pdf with
% the best fit page size

clear axisLimits cPlotCount F figureSize fileList filename filenamePrefix ...
    fileNo folder frameNo graphSelection graphTitle graphTitleList h ...
    makeMovie matShape movieName numTimeStepLimit numTimeStepsToGraph   ...
    plotFreq quiverScaleOnX removeBoundaryVal res singleFrame ...
    setViewingAngle numFrames solnNo t titleLine2 tt v viewingAngle ...
    zAxisLimits
