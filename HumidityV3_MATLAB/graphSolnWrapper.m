function graphSolnWrapper(projectPath, solnFileFullPath, param, solnName, plotName, stepNo, ...
    graphGhostCells, graphContourPlots, viewAngle, saveToPDF, contourProportion, colmapName)
%GRAPHSOLNWRAPPER  Wrapper function for graphing numerical solutions, numerical errors, or exact solutions

saveSurfPlotToPDF = false;

% Parsing parameters
if nargin < 7
    graphGhostCells = false;
end
if nargin < 8
    graphContourPlots = false;
end
if nargin < 9
    viewAngle = false;
end
if nargin < 10
    saveToPDF = false;
end
if nargin < 11
    contourProportion = 0.5;
end
if nargin < 12
    colmapName = 'default';  % or 'gray'
end

if ~exist(solnFileFullPath, 'file')
    fprintf("File does not exist (yet)!\n");
    return
end
fprintf("Opening file: " + solnFileFullPath + "\n");

% Mesh grid: cell centers
centersX = csvread(projectPath + "Output/CellCenters_X.csv");
centersP = csvread(projectPath + "Output/CellCenters_P.csv");

% Generate graph titles
titleLine1 = "Plot of " + plotName + " " + solnName;
plotTime = param.Dt * stepNo;
displayPlotTime = plotTime;
titleLine2 = "Time = " + num2str(displayPlotTime) + "s";

if graphGhostCells
    % Graph all solution values, with ghost cell values
    if viewAngle == false
        fig = graphSoln(solnFileFullPath, centersX, centersP, titleLine1, titleLine2);
    else
        fig = graphSoln(solnFileFullPath, centersX, centersP, titleLine1, titleLine2, viewAngle);
    end
else
    % No ghost cells
    centersX_noGhost = centersX(2:end-1, 2:end-1);
    centersP_noGhost = centersP(2:end-1, 2:end-1);
    if viewAngle == false
        fig = graphSoln(solnFileFullPath, centersX_noGhost, centersP_noGhost, titleLine1, titleLine2);
    else
        fig = graphSoln(solnFileFullPath, centersX_noGhost, centersP_noGhost, titleLine1, titleLine2, viewAngle);
    end
end

% Print pdf
if saveToPDF && saveSurfPlotToPDF
    filename = "Output/" + solnName + "_" + num2str(displayPlotTime) + "s";
    saveas(fig, char(filename), 'pdf');
end

% Graph contour plots
if graphContourPlots
    contourLevels = selectContourLevels(solnName);
    figure
    if graphGhostCells
        fig = graphContour(solnFileFullPath, centersX, centersP, contourLevels, titleLine1, titleLine2, ...
            viewAngle, contourProportion, colmapName);
    else
        fig = graphContour(solnFileFullPath, centersX_noGhost, centersP_noGhost, contourLevels, titleLine1, titleLine2, ...
            viewAngle, contourProportion, colmapName);
    end
    if saveToPDF
        filename = "Output/" + solnName + "_" + num2str(displayPlotTime) + "s_contour";
        saveas(fig, char(filename), 'pdf');
    end
end

end

