%% Examples: This file provides simple examples of surface plots of the numerical results
%  Following examples include plotting of numerical solutions, numerical errors, and exact solutions

%% Example 1: Graph solutions
clear; clc

% ===== ===== ===== =====
% Settings
% ===== ===== ===== =====

projectPath = "/home/chuck/git/Numerical-Repo/HumidityV3/";  % Path to the outmost folder
resultFolder = "MovieFrames/";  % Commonly used: "MovieFrames/" or "Output/"

plotName = "Numerical Solution";  % "Solution", "Error", or "Exact Solution"
solnName = "T";
stepNo = 1001;

graphGhostCells = false;

% ===== ===== ===== =====
% Graph solution
% ===== ===== ===== =====

param = readParam(projectPath + "Output/Param.csv");
solnFileFullPath = genSolnFileFullPath(projectPath, resultFolder, solnName, stepNo);
graphSolnWrapper(projectPath, solnFileFullPath, param, solnName, plotName, stepNo, graphGhostCells)


%% Example 2: Graph errors and exact solutions

clear; clc

% ===== ===== ===== =====
% Settings
% ===== ===== ===== =====

% Path to the outmost folder
projectPath = "/home/chuck/git/Numerical-Repo/HumidityV3/";
resultFolder = "Output/";  % Commonly used: "MovieFrames/" or "Output/"

plotName = "Exact Solution";  % "Solution", "Error", or "Exact Solution"
solnName = "T";

graphGhostCells = false;

% ===== ===== ===== =====
% Graph Errors
% ===== ===== ===== =====

param = readParam(projectPath + "Output/Param.csv");
stepNo = param.numTimeStep;
if plotName == "Error"
    solnFileFullPath = projectPath + resultFolder + solnName + "_err.csv";
elseif plotName == "Exact Solution"
    solnFileFullPath = projectPath + resultFolder + solnName + "_exact.csv";
end
graphSolnWrapper(projectPath, solnFileFullPath, param, solnName, plotName, stepNo, graphGhostCells)


%% Example 3: Graph contour plots of solutions

clearAllScp

% ===== ===== ===== =====
% Settings
% ===== ===== ===== =====

projectPath = "~/git/Numerical-Repo/HumidityV3/";  % Path to the outmost folder
resultFolder = "MovieFrames/";  % Commonly used: "MovieFrames/" or "Output/"

plotName = "Numerical Solution";  % "Solution", "Error", or "Exact Solution"
solnName = "q";
stepNo = 19500;

saveToPDF = false;
viewAngle = [0, -90];
graphGhostCells = false;
graphContourPlots = true;

% ===== ===== ===== =====
% Graph solution
% ===== ===== ===== =====

param = readParam(projectPath + "Output/Param.csv");
solnFileFullPath = genSolnFileFullPath(projectPath, resultFolder, solnName, stepNo);
graphSolnWrapper(projectPath, solnFileFullPath, param, solnName, plotName, stepNo, ...
    graphGhostCells, graphContourPlots, viewAngle, saveToPDF)


%% Wrapper function for graphing numerical solutions, numerical errors, or exact solutions

function graphSolnWrapper(projectPath, solnFileFullPath, param, solnName, plotName, stepNo, ...
    graphGhostCells, graphContourPlots, viewAngle, saveToPDF)

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

contourProportion = 1;

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
titleLine2 = "Time = " + num2str(plotTime) + "s";

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
if saveToPDF
    filename = "Output/" + solnName + "_" + num2str(plotTime) + "s";
    saveas(fig, char(filename), 'pdf');
end

% Graph contour plots
if graphContourPlots
    contourLevels = selectContourLevels(solnName);
    figure
    if graphGhostCells
        fig = graphContour(solnFileFullPath, centersX, centersP, contourLevels, titleLine1, titleLine2, ...
            viewAngle, contourProportion);
    else
        fig = graphContour(solnFileFullPath, centersX_noGhost, centersP_noGhost, contourLevels, titleLine1, titleLine2, ...
            viewAngle, contourProportion);
    end
    if saveToPDF
        filename = "Output/" + solnName + "_" + num2str(plotTime) + "s_contour";
        saveas(fig, char(filename), 'pdf');
    end
end

end






















