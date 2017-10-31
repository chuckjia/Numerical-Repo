% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== 
% Read Grid Points From File
% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== %

% To run this file, it is required to run getPar_sct first.

folder = "Results/";
matDim = [numCellsP + 1, numCellsX + 1];  % reshape is column major
meshGridX = reshape(getVecFromFile_fcn(folder, "meshGridX"), matDim);  
meshGridP = reshape(getVecFromFile_fcn(folder, "meshGridP"), matDim);

clear folder matDim

