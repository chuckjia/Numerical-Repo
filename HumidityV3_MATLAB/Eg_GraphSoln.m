% This file is an example of graphing surface and contour plots of numerical solutions.
%

clearAllScp

% ===== ===== ===== =====
% Settings
% ===== ===== ===== =====

solnName = "q";
stepNo = -1;  % A value of -1 indicates graph the latest solution

projectPath = "~/git/Numerical-Repo/HumidityV3/";  % Path to the outmost folder
resultFolder = "MovieFrames/";  % Commonly used: "MovieFrames/" or "Output/"
plotName = "Numerical Solution";  % "Solution", "Error", or "Exact Solution"

saveToPDF = false;  % Whether to save result to PDF files
viewAngle = [0, -90];  % Viewing angle
graphGhostCells = false;  % Whether to graph ghost cells
graphContourPlots = true;  % Whether to graph contour plots

% Setting some of the values: Do Not Change
stepNo = genActualStepNo(stepNo, projectPath + resultFolder, solnName);

% ===== ===== ===== =====
% Graph solution
% ===== ===== ===== =====

param = readParam(projectPath + "Output/Param.csv");
solnFileFullPath = genSolnFileFullPath(projectPath, resultFolder, solnName, stepNo);
graphSolnWrapper(projectPath, solnFileFullPath, param, solnName, plotName, stepNo, ...
    graphGhostCells, graphContourPlots, viewAngle, saveToPDF)









