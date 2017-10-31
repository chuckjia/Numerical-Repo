% 
% Test_QuadCell_coef.m
%
% This file contains the tests on the calculation of interpolation
% coefficients for the quadrilateral cells. More specifically, these tests  
% are about the linear solver in the c code for (3.13).
%

clear; clc

% Read data from files
getPar_sct;
getCellCenters_sct;
cellCentersX = cellCentersX.'; cellCentersP = cellCentersP.';
getMeshGrids_sct;
meshGridX = meshGridX.'; meshGridP = meshGridP.';

% Read interpolation coefficients from file
folder = "Results/";
matShape = [Np + 1, Nx + 1];
a2Mat = reshape(getVecFromFile_fcn(folder, "a2_quadCell"), matShape).';
a3Mat = reshape(getVecFromFile_fcn(folder, "a3_quadCell"), matShape).';
a4Mat = reshape(getVecFromFile_fcn(folder, "a4_quadCell"), matShape).';
resFromC = zeros(Nx + 1, Np + 1, 3);
resFromC(:,:,1) = a2Mat;
resFromC(:,:,2) = a3Mat;
resFromC(:,:,3) = a4Mat;

a1 = 0.25;
res = zeros(Nx + 1, Np + 1, 3);
err = zeros(Nx + 1, Np + 1);
for i = 0:Nx
    for j = 1:Np
        i_ind = i + 1; j_ind = j + 1;
        coefMat = [...
            cellCentersX(i_ind + 1, j_ind), cellCentersX(i_ind, j_ind + 1), cellCentersX(i_ind + 1, j_ind + 1); ...
            cellCentersP(i_ind + 1, j_ind), cellCentersP(i_ind, j_ind + 1), cellCentersP(i_ind + 1, j_ind + 1); ...
            1, 1, 1];
        bMat = [ ...
            meshGridX(i_ind + 1, j_ind + 1) - a1 * cellCentersX(i_ind, j_ind); ...
            meshGridP(i_ind + 1, j_ind + 1) - a1 * cellCentersP(i_ind, j_ind); ...
            1 - a1];
        resThis = coefMat \ bMat;
        res(i_ind, j_ind, :) = resThis;
        resFromCThis = reshape(resFromC(i_ind, j_ind, :), [3, 1]);
        err(i_ind, j_ind) = norm(resThis - resFromCThis, 2);
    end
end

reshapeDim = [1, (Nx + 1) * (Np + 1)];
maxErr = max(reshape(err, reshapeDim));
fprintf("The largest error is of norm %1.8e\n", maxErr);

clear folder reshapeDim maxErr
