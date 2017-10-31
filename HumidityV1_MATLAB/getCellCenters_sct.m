% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
% Read Mesh: Cell Centers
% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== %

% x-coordinates and p-coordinates for the cell centers

% To run this file, getPar_sct is required to be executed first

% Note that reshape() method works in column-major order. Thus the shape is 
% [numCellsP, numCellsX], the transpose from the matrix in c

matShape = [numCellsP, numCellsX];
folder = "Results/";
cellCentersX = reshape(getVecFromFile_fcn(folder, "cellCentersX"), matShape);
cellCentersP = reshape(getVecFromFile_fcn(folder, "cellCentersP"), matShape);

clear matShape folder

