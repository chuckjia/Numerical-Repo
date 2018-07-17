%% Examples: This file provides simple examples of surface plots of the numerical results
%  Following examples include plotting of numerical solutions, numerical errors, and exact solutions

%% Example 1: Graph solutions

clearAllScp

% ===== ===== ===== =====
% Settings
% ===== ===== ===== =====

projectPath = "/home/chuck/git/Numerical-Repo/HumidityV3/";  % Path to the outmost folder
resultFolder = "MovieFrames/";  % Commonly used: "MovieFrames/" or "Output/"

plotName = "Numerical Solution";  % "Solution", "Error", or "Exact Solution"
solnName = "T";
stepNo = 1000;

graphGhostCells = false;

% ===== ===== ===== =====
% Graph solution
% ===== ===== ===== =====

param = readParam(projectPath + "Output/Param.csv");
solnFileFullPath = genSolnFileFullPath(projectPath, resultFolder, solnName, stepNo);
graphSolnWrapper(projectPath, solnFileFullPath, param, solnName, plotName, stepNo, graphGhostCells)


%% Example 2: Graph errors and exact solutions

clearAllScp

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

