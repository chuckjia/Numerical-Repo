% Graph_PhysSim.m - Graph numerical solutions for the physical simulation
%
% Author: Chuck Jia
% Created on: Dec 20, 2017

clear; clc

for numThou = 0:10
clearvars -except numThou

folder = "Sim10/"

path = "/Users/chuckjia/Documents/Workspace/DataStorage/Humidity/";
folder_getParScp = path + folder + "Norm/"; getPar_scp; getCellCenters_scp;

% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== %
% Plot Settings
% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== %

solnName = 'w';
% numThou = 20;
timeToPlot = 1000 * numThou

zCoordPlot = false;
plotBoundary = true;
plotSurf = true;
plotContourf = true;

sourceFolder = path + folder + "Snapshots/";

saveFile = true;

saveAsPDF = true;
saveAsFig = true;
if (saveFile == false)
    saveAsPDF = false;
    saveAsFig = false;
end

% ----- ----- ----- ----- ----- %
% Process paramters
% ----- ----- ----- ----- ----- %

switch solnName
    case "T"
        % contourLevels = [267, 273, 279, 284, 290, 293, 296, 299, 350];  % For height 250
        % contourLevels = [270, 284, 290, 293, 295, 296, 298, 299, 350];  % For height 250
        % contourLevels = [267, 273, 279, 284, 290, 293, 295, 297];
         contourLevels = [190, 267, 273, 279, 284, 290, 293, 295, 297];
    case "q"
        % contourLevels = 0.01 * [0.6, 0.8, 0.9, 1.1, 1.2, 1.3, 1.4, 1.5, 1.59, 1.7];  % p coord, height = 250
        % contourLevels = 1e-3 * [1, 3, 5, 8.83, 10.7, 12.5, 13.5, 14.5, 15.5];  % z coord, height = 250
         contourLevels = 1e-3 * [1, 3, 5, 8.83, 10.7, 12.5, 13.5, 14.5];
    case "u"
        % contourLevels = [5, 5.6, 5.7, 5.8, 6, 6.5, 7, 8, 9, 9.7, 10];
        contourLevels = [2, 3, 4, 6, 6.8, 7.5, 8, 9, 10, 11];
    case "w"
        % contourLevels = [-1.4, -0.4, -0.3, -0.2, -0.1, -0.05, 0, 0.05, 0.1, 0.2, 0.3, 0.4];
        contourLevels = [-1.4, -0.8, -0.4, -0.1, 0, 0.1, 0.4, 0.6];
    otherwise
        disp("Not a valid solnName!\n")
end

% ----- ----- ----- ----- ----- %
% View parameters
% ----- ----- ----- ----- ----- %

if (zCoordPlot)
    yCoordName = 'z coordinate';
    yAxisRange = [0, 6000]; nearMountPortion = 12.5 / 25;
    % yAxisRange = [0, 10000]; nearMountPortion = 1 / numCellsP;
    % yAxisRange = [0, convPtoZ_fcn(pA)]; nearMountPortion = 1 / numCellsP;
else  % p Coordinates
    yCoordName = 'p coordinate';
    yAxisRange = [pA, pB]; nearMountPortion = 1 / numCellsP;
    % yAxisRange = [600, pB]; nearMountPortion = 20 / 25;  % yAxisRange from 650 in original simulation
end

% viewAngle = [0, -90];
viewAngle = [0, 90];
nearBoundCells = 0;

% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== %
% Graph Selected Plots
% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== %

% ----- ----- ----- ----- ----- %
% Process data from C program
% ----- ----- ----- ----- ----- %

matShape = [numCellsX, numCellsP];
% Choose whether to plot the boundary or not
if (~plotBoundary)
    cellCentersX = rmBDVal_fcn(cellCentersX);
    cellCentersP = rmBDVal_fcn(cellCentersP);
end
if (zCoordPlot)
    cellCentersP = convPtoZ_fcn(cellCentersP);
end

fprintf("Plotting\n");

% Read numerical solutions/errors from file
filename = solnName + "_soln_" + int2str(floor(timeToPlot / Dt));
titleText = solnName + " at t = " + timeToPlot + "s";
% Reshape matrix read from file
res = reshape(getVecFromFile_fcn(sourceFolder, filename), matShape);
if (~plotBoundary) res = rmBDVal_fcn(res); end
% Plot the solution and error
pFloor = floor(numCellsP * nearMountPortion);
pCeiling = size(cellCentersP, 1) - nearBoundCells;
x_coords = cellCentersX(pFloor:pCeiling,:);
p_coords = cellCentersP(pFloor:pCeiling,:);
res = res(pFloor:pCeiling,:);
if (plotSurf)
    figure
    surf(x_coords, p_coords, res);
    xlim([x0, xf]); ylim(yAxisRange); view(viewAngle); colormap(jet)
    title(titleText); xlabel('x coordinates'); ylabel(yCoordName);
    
    filename_print = strcat('Results/Surf_', solnName, '_at_', int2str(timeToPlot), 's');
    if (saveAsPDF) print(filename_print, '-dpdf', '-bestfit'); end
    if (saveAsFig) savefig(filename_print); end
end
if (plotContourf)
    figure;
    contourf(x_coords, p_coords, res, contourLevels);
    xlim([x0, xf]); ylim(yAxisRange); view(viewAngle); colormap(jet)
    title(titleText); xlabel('x coordinates'); ylabel(yCoordName);
    colorbar('eastoutside')
    filename_print = strcat('Results/Contour_', solnName, '_at_', int2str(timeToPlot), 's');
    if (saveAsPDF) print(filename_print, '-dpdf', '-bestfit'); end
    if (saveAsFig) savefig(filename_print); end
end

clear titleLine2 matShape fileList graphTitleList graphList

end
