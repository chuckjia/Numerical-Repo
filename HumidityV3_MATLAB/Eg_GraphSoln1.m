clear; clc

% Path to the outmost folder
path = "/Users/chuckjia/Documents/Workspace/Git/Numerical-Repo/HumidityV3/";
solnFilename = path + "MovieFrames/T_1000.csv";

% Mesh grid: cell centers
centersX = csvread(path + "Output/CellCenters_X.csv");
centersP = csvread(path + "Output/CellCenters_P.csv");

% Generate graph titles
titleLine1 = "Plot of Solution u";
titleLine2 = "Time = 1s";

% Graph all solution values, with ghost cell values
graphSoln(solnFilename, centersX, centersP, titleLine1, titleLine2);

% % No ghost cells
% centersX_noGhost = centersX(2:end-1, 2:end-1);
% centersP_noGhost = centersP(2:end-1, 2:end-1);
% graphSoln(solnFilename, centersX_noGhost, centersP_noGhost, titleLine1, titleLine2);