% This file is an example of graphing surface and contour plots of numerical solutions.
%

clearAllScp
projectPath = defaultProjectPath();  % Path to the outmost folder
projectPath = "~/Documents/Workspace/DataStorage/Humidity/HumidityV3/2018_11_29_15_39_25_good";
projectPath = genFolderPathName(projectPath);

% ===== ===== ===== =====
% Settings
% ===== ===== ===== ===== 

solnName = "q";
stepNo = -1;  % A value of -1 indicates graph the latest solution

resultFolder = "MovieFrames/";  % Commonly used: "MovieFrames/" or "Output/"
plotName = "Numerical Solution";  % "Solution", "Error", or "Exact Solution"

saveToPDF = true;  % Whether to save result to PDF files
viewAngle = [0, -90];  % Viewing angle
graphGhostCells = false;  % Whether to graph ghost cells
graphContourPlots = true;  % Whether to graph contour plots

contourProportion = 0.5;
myGrayMap = linspace(0.05, 0.95, 15);  myGrayMap = repmat(myGrayMap', 1, 3);
colmapName = 'default'; % myGrayMap, 'default' or 'gray'

% Setting some of the values: Do Not Change
stepNo = genActualStepNo(stepNo, projectPath + resultFolder, solnName);

% ===== ===== ===== =====
% Graph solution
% ===== ===== ===== =====

param = readParam(projectPath + "Output/Param.csv");
solnFileFullPath = genSolnFileFullPath(projectPath, resultFolder, solnName, stepNo);
graphSolnWrapper(projectPath, solnFileFullPath, param, solnName, plotName, stepNo, ...
    graphGhostCells, graphContourPlots, viewAngle, saveToPDF, contourProportion, colmapName)


