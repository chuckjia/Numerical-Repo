clear; clc

% Settings
solnName = "T";

% Path and solution file names
path = "/Users/chuckjia/Documents/Workspace/Git/Numerical-Repo/HumidityV3/";  % Path to the outmost folder

% Mesh and paramters
param = readParam(path + "Output/Param.csv"); 
centersX = csvread(path + "Output/CellCenters_X.csv"); centersX = centersX(2:end-1, 2:end-1);
centersP = csvread(path + "Output/CellCenters_P.csv"); centersP = centersP(2:end-1, 2:end-1);

steps = 100 * (1:200);
F = genMovie(solnName, path + "MovieFrames/", param, centersX, centersP, steps, true);
% slowMotion(F, 10, char("Output/q_slow.avi"));


