% Graph_PhysSimuln.m - Graph numerical solutions for the physical simulation
%
% Author: Chuck Jia
% Created on: Dec 20, 2017

clear; clc

% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== %
% Plot Settings
% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== %

solnName = 'T';
timeToPlot = 1000 * 10;
plotBoundary = false;
plotSurf = false;
plotContourf = true;
source_folder = "Storage/NewSim/";

switch solnName
    case "T"
        contourLevels = [280, 284, 290, 294, 296, 298, 299, 300, 308];
    case "q"
        contourLevels = 0.01 * ...
            [0.070, 0.883, 1.07, 1.25, 1.35, 1.45, 1.55, 1.6, 1.7, 1.8] ;
    case "u"
        % contourLevels = [5, 5.6, 5.7, 5.8, 6, 6.5, 7, 8, 9, 9.7, 10];
        contourLevels = [4, 6, 7.5, 8, 9, 10, 11];
    case "w"
        contourLevels = [-0.4, -0.3, -0.2, -0.1, -0.05, 0, ...  
            0.05, 0.1, 0.2, 0.3, 0.4];
    otherwise
        disp("Not a valid solnName!\n")
end

% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== %
% Graph Selected Plots
% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== %

getPar_scp;  % Read and calculate parameters
getCellCenters_scp;  % Read cell centers
% y_axis_range = [pA, pB]; near_mount_portion = 1 / numCellsP;
y_axis_range = [600, pB]; near_mount_portion = 20 / 25;  % y_axis_range from 650 in original simulation
viewAngle = [0, -90];

matShape = [numCellsX, numCellsP];
% Choose whether to plot the boundary or not
if (~plotBoundary)
    cellCentersX = rmBDVal_fcn(cellCentersX);
    cellCentersP = rmBDVal_fcn(cellCentersP);
end
fprintf("Plotting\n");

% Read numerical solutions/errors from file
filename = solnName + "_soln_" + int2str(floor(timeToPlot / Dt));
titleText = solnName + " at t = " + timeToPlot + "s";

res = reshape(getVecFromFile_fcn(source_folder, filename), matShape);
if (~plotBoundary)
    res = rmBDVal_fcn(res);
end
% Plot the solution and error
p_floor = floor(numCellsP * near_mount_portion); 
p_ceiling = numCellsP - 2 - 3;
x_coords = cellCentersX(p_floor:p_ceiling,:);
p_coords = cellCentersP(p_floor:p_ceiling,:);
res = res(p_floor:p_ceiling,:);
if (plotSurf)
    figure
    surf(x_coords, p_coords, res);
    xlim([x0, xf]); ylim(y_axis_range);
    view(viewAngle)
    title(titleText)
    xlabel('x coordinates'); ylabel('p coordinates');
end
if (plotContourf)
    fig = figure;
    contourf(x_coords, p_coords, res, contourLevels);
    xlim([x0, xf]); ylim(y_axis_range);
    % colormap(hot)
    view(viewAngle)
    title(titleText)
    colorbar('eastoutside')
    % Add titles and labels
    xlabel('x coordinates'); ylabel('p coordinates');
    filename_print = strcat('Results/', solnName, '_at_', int2str(timeToPlot), 's')
    print(filename_print, '-dpdf', '-bestfit')
end


clear titleLine2 matShape fileList graphTitleList graphList
