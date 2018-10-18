% This file is an example of graphing surface and contour plots of numerical solutions.
%

clearAllScp
% projectPath = "~/Documents/Workspace/git/Numerical-Repo/HumidityV3/";  % Path to the outmost folder
% projectPath = "~/git/Numerical-Repo/HumidityV3/";  % Path on Ubuntu
projectPath = "~/Documents/Workspace/DataStorage/Humidity/HumidityV3/2018_10_17_11_1_7";

projectPath = fixProjectPathName(projectPath);

% ===== ===== ===== =====
% Settings
% ===== ===== ===== =====

solnName = "q";
stepNo = 40000;  % A value of -1 indicates graph the latest solution

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


%% 

function projectPath = fixProjectPathName(projectPath)
    if ischar(projectPath)
        lastChar = projectPath(end);
        projectPath = string(projectPath);
    else
        projectPath_temp = char(projectPath);
        lastChar = projectPath_temp(end);
    end
    
    if lastChar ~= '/'
        projectPath = projectPath + "/";
    end
end






