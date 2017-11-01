% Graph_Mesh.m - Plot the mesh grid and cell centers
% 
% This script plots the mesh grid and cell centers from the humidity
% calculations.
%
% Author: Chuck Jia
% Created on: Oct 20, 2017

clear; clc

% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== %
% Plot Settings
% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== %

% Control whether to plot cell centers
plotCenter = true;

% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== %
% Plot The Mesh And Cell Centers
% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== %

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
