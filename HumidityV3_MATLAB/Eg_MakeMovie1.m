% This file is an example of making a movie of the numerical simulation.
%

clearAllScp

% ===== ===== ===== ===== 
% Settings
% ===== ===== ===== ===== 

solnName = "q";
steps = 0:500:40000;  % Vector of all the step numbers to be included in movie
% Path and solution file names
% projectPath = "~/Documents/Workspace/Git/Numerical-Repo/HumidityV3/";  % Path to the outermost folder
projectPath = "~/Documents/Workspace/DataStorage/Humidity/HumidityV3/2018_11_29_15_39_25_good/";
outputFilename = "Output/" + solnName + "_slow.avi";

savefile = false;

% ===== ===== ===== ===== 
% Generate Movies
% ===== ===== ===== ===== 

% Mesh and paramters
param = readParam(projectPath + "Output/Param.csv"); 
centersX = csvread(projectPath + "Output/CellCenters_X.csv"); centersX = centersX(2:end-1, 2:end-1);
centersP = csvread(projectPath + "Output/CellCenters_P.csv"); centersP = centersP(2:end-1, 2:end-1);

F = genMovie(solnName, projectPath + "MovieFrames/", param, centersX, centersP, steps, false);
fprintf("Parameters used: mesh size = %dx%d, averaging frequency = every %d steps.\n", param.Nx, param.Np, param.aveFreq);
if savefile
    slowMotion(F, 5, char(outputFilename));  
end


