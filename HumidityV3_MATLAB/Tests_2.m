
%% 

clear; clc; close all

projectPath = "~/Documents/Workspace/Git/Numerical-Repo/HumidityV3/";  % Path to the outermost folder

% Mesh and paramters
param_0 = readParam(projectPath + "Output/Param.csv"); 
centersX_0 = csvread(projectPath + "Output/CellCenters_X.csv"); centersX_0 = centersX_0(2:end-1, 2:end-1);
centersP_0 = csvread(projectPath + "Output/CellCenters_P.csv"); centersP_0 = centersP_0(2:end-1, 2:end-1);
meshGridX_0 = csvread(projectPath + "Output/MeshGrid_X.csv");
meshGridP_0 = csvread(projectPath + "Output/MeshGrid_P.csv");

projectPath = "~/Documents/Workspace/Eclipse/JIASpace/HumidityTest/";

% Mesh and paramters
param_1 = readParam(projectPath + "Output/Param.csv"); 
centersX_1 = csvread(projectPath + "Output/CellCenters_X.csv"); centersX_1 = centersX_1(2:end-1, 2:end-1);
centersP_1 = csvread(projectPath + "Output/CellCenters_P.csv"); centersP_1 = centersP_1(2:end-1, 2:end-1);
meshGridX_1 = csvread(projectPath + "Output/MeshGrid_X.csv");
meshGridP_1 = csvread(projectPath + "Output/MeshGrid_P.csv");

% projectPath = "~/Documents/Workspace/Eclipse/JIASpace/HumidityTempTest/";
% % Mesh and paramters
% % param_1 = readParam(projectPath + "Output/Param.csv"); 
% centersX_2 = csvread(projectPath + "Output/CellCenters_X.csv"); centersX_2 = centersX_2(2:end-1, 2:end-1);
% centersP_2 = csvread(projectPath + "Output/CellCenters_P.csv"); centersP_2 = centersP_2(2:end-1, 2:end-1);
% meshGridX_2 = csvread(projectPath + "Output/MeshGrid_X.csv");
% meshGridP_2 = csvread(projectPath + "Output/MeshGrid_P.csv");

max(reshape(centersX_1 - centersX_0, [], 1))
max(reshape(centersP_1 - centersP_0, [], 1))

max(reshape(meshGridX_1 - meshGridX_0, [], 1))
max(reshape(meshGridP_1 - meshGridP_0, [], 1))

%% Compare values of 4 solutions

clear; clc; close all

solnName = 'T';
stepNo = 0;

projectPath = "~/Documents/Workspace/Git/Numerical-Repo/HumidityV3/";  % Path to the outermost folder
filename = projectPath + "MovieFrames/" + solnName + "_" + num2str(stepNo) + ".csv";
sl_0 = csvread(filename);  

projectPath = "~/Documents/Workspace/Eclipse/JIASpace/HumidityTest/";
filename = projectPath + "MovieFrames/" + solnName + "_" + num2str(stepNo) + ".csv";
sl_1 = csvread(filename);

abs(sl_0 - sl_1)


%% Compare gradh values

clear; clc; close all

projectPath = "~/Documents/Workspace/Git/Numerical-Repo/HumidityV3/";  % Path to the outermost folder
gradhU_x_0 = csvread(projectPath + "Output/gradhU_x.csv"); 
% quad_a2_0 = csvread(projectPath + "Output/a2_quad.csv");
% quad_a3_0 = csvread(projectPath + "Output/a3_quad.csv");
% quad_a4_0 = csvread(projectPath + "Output/a4_quad.csv");
% quad_MInv12_0 = csvread(projectPath + "Output/e12_MInv_quad.csv");
getCellTopRightU_0 = csvread(projectPath + "Output/CellTopRightU.csv");

projectPath = "~/Documents/Workspace/Eclipse/JIASpace/HumidityTest/";
filename = projectPath + "Output/gradhU_x.csv";
gradhU_x_1 = csvread(filename);
% quad_a2_1 = csvread(projectPath + "Output/a2_quad.csv");
% quad_a3_1 = csvread(projectPath + "Output/a3_quad.csv");
% quad_a4_1 = csvread(projectPath + "Output/a4_quad.csv");
% quad_MInv12_1 = csvread(projectPath + "Output/e12_MInv_quad.csv");
getCellTopRightU_1 = csvread(projectPath + "Output/CellTopRightU.csv");

% gradhU_x_0 - gradhU_x_1
% quad_a3_0 - quad_a3_1
% quad_a4_0 - quad_a4_1
% max(max(abs(quad_MInv12_0 - quad_MInv12_1)))
getCellTopRightU_0 - getCellTopRightU_1




%% Compare fluxes

clear; clc; close all

fluxName = 'GG_T';

projectPath = "~/Documents/Workspace/Git/Numerical-Repo/HumidityV3/";  % Path to the outermost folder
filename = projectPath + "Output/" + fluxName + ".csv";
flux_0 = csvread(filename);  

projectPath = "~/Documents/Workspace/Eclipse/JIASpace/HumidityTest/";
filename = projectPath + "Output/" + fluxName + ".csv";
flux_1 = csvread(filename);

abs(flux_0 - flux_1)


%%

%% Compare values of 4 solutions

clear; clc; close all

solnName = 'w';
stepNo = 0;

projectPath = "~/Documents/Workspace/Eclipse/JIASpace/HumV1Broken/";  % Path to the outermost folder
filename = projectPath + "MovieFrames/" + solnName + "_" + num2str(stepNo) + ".csv";
sl_0 = csvread(filename);  

projectPath = "~/Documents/Workspace/Eclipse/JIASpace/HumV1Good/";
filename = projectPath + "MovieFrames/" + solnName + "_" + num2str(stepNo) + ".csv";
sl_1 = csvread(filename);

max(abs(sl_0 - sl_1))











