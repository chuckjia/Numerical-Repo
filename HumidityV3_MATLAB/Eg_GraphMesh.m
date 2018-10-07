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

% hold on 
% x0 = 0;  xf = 75000;  Nx = 10;  Np = 10;  Dx = (xf - x0) / Nx;
% pB = @(x) 1000 - 250 * exp(-((x - 37500) / 6000) .^ 2);
% npt = 200;
% % x_vec = linspace(x0-Dx, xf+Dx, npt);
% x_vec = linspace(x0, xf, npt);
% plot3(x_vec, pB(x_vec), zeros(1, npt),'-', 'Color', 0.5*[1,1,1], 'LineWidth', 2);
% hold off