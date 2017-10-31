clear; clc

% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
% Test Matrices in the Quadrilateral Cells
% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== %

% This file contains the tests the calculation of M_{i,j+1/2} in (3.15)

% Read data from files
getPar_sct;
getCellCenters_sct;
cellCentersX = cellCentersX.'; cellCentersP = cellCentersP.';
getMeshGrids_sct;
meshGridX = meshGridX.'; meshGridP = meshGridP.';

% Read matrix elements from file
folder = "Results/";
matShape = [Np + 1, Nx + 1];
e12_mat = reshape(getVecFromFile_fcn(folder, "e12_diagMat_quadCell"), matShape).';
e22_mat = reshape(getVecFromFile_fcn(folder, "e22_diagMat_quadCell"), matShape).';

e11_val = getVecFromFile_fcn(folder, "e11e21_diagMat_quadCell");
e21_val = e11_val(2); e11_val = e11_val(1);

err = ones((Nx + 1) * Np, 1);
count = 0;
for i = 0:Nx
    for j = 1:Np
        count = count + 1;
        i_m = i + 1; j_m = j + 1;
        M = [
            meshGridX(i_m + 1, j_m + 1) - meshGridX(i_m, j_m + 1), ...
            meshGridP(i_m + 1, j_m + 1) - meshGridP(i_m, j_m + 1);
            cellCentersX(i_m, j_m + 1) - cellCentersX(i_m, j_m), ...
            cellCentersP(i_m, j_m + 1) - cellCentersP(i_m, j_m)
            ];
        M_inv = inv(M);
        M_fromC = [e11_val, e12_mat(i_m, j_m); e21_val, e22_mat(i_m, j_m)];
        err(count) = norm(reshape(M_inv - M_fromC, [1, 4]), 2);
    end
end

fprintf("The max error has norm %1.10e\n", max(err));








