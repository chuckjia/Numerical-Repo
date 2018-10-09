clear; clc; close all
tic

% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== 
% Set Model Parameters
% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== 

modelNo = 0;
param = setModelParam();
[x0, xf, pA, Nx, Np, Dt, Nt, t_end, AveRate, MovieFrameRate, NumMsg] = param.extractParam();
numCellsX = Nx + 2; numCellsP = Np + 2;
Dx = (xf - x0) / Nx;

% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== 
% Build Mesh
% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== 

[meshGridX, meshGridP, DpVec] = buildMesh(param);
[cellCentersX, cellCentersP] = buildCellCenters(meshGridX, meshGridP);

a1 = 0.25;
[a2mat, a3mat, a4mat] = buildQuadCells(meshGridX, meshGridP, cellCentersX, cellCentersP, a1);
[M_e12, M_e22] = buildQuadM(meshGridP, cellCentersP);

% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== 
% Initialize Solutions
% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====

[T_, q_, u_] = enforceInitCond(param, cellCentersX, cellCentersP, modelNo);
w_ = zeros(numCellsX, numCellsP);
phix_ = zeros(numCellsX, numCellsP);
gradxU = calcGradxU(Dx, u_, a1, a2mat, a3mat, a4mat, M_e12, M_e22);



toc











