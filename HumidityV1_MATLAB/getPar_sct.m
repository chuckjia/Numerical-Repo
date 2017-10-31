% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
% Read Parameters From File 
% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== %

folder = "Results/";
mat = getVecFromFile_fcn(folder, "par");
x0 = mat(1); xf = mat(2);
pA = mat(3); pBx0 = mat(4); pBxf = mat(5); pB = pBx0;
Nx = mat(6); Np = mat(7);
Dt = mat(8); numTimeSteps = mat(9);

% Calculate related parameters
Dx = (xf - x0) / Nx; Dp = (pB - pA) / Np;
numCellsX = Nx + 2; numCellsP = Np + 2;

clear mat folder