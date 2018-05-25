clear; clc
folder_path = "/Users/chuckjia/Documents/Workspace/Git/Numerical-Repo/HumidityV3/Output/";
gridFile_X = folder_path + "MeshGrid_X.csv";
gridFile_P = folder_path + "MeshGrid_P.csv";
centersFile_X = folder_path + "CellCenters_X.csv";
centersFile_P = folder_path + "CellCenters_P.csv";

% Graph mesh without marking cell centers
% figure;
% graphMesh(grid_x_file, grid_p_file);

% Graph mesh with cell centers
% figure;
 graphMesh(gridFile_X, gridFile_P, centersFile_X, centersFile_P);