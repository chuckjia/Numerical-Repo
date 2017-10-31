% 
% Plot The Mesh And Cell Centers
%

clear; clc

% Control plotting cell centers
plotCenter = true;

% Read data from file
getPar_sct;  % Read and calculate parameters
getMeshGrids_sct;  % Read mesh grids from file

% Plot the mesh grid
mesh(meshGridX, meshGridP, zeros(numCellsX + 1, numCellsP + 1))
view([0, 0, 1])
hold on

% Plot Cell Centers
if (plotCenter)
    getCellCenters_sct;
    plot3(cellCentersX, cellCentersP, zeros(numCellsP, numCellsX), ...
        '.', 'MarkerSize', 10)
end
hold off
