% getPar_scp.m - Read parameters from file
%
% Read parameter from file for later use
%

if (exist(folder_getParScp) == 0)
    folder_getParScp = "/Users/chuckjia/Documents/Workspace/Git/Numerical-Repo/HumidityV1/Results/";
end
mat_getParScp = getVecFromFile_fcn(folder_getParScp, "par");
x0 = mat_getParScp(1); xf = mat_getParScp(2);
pA = mat_getParScp(3); pBx0 = mat_getParScp(4); pBxf = mat_getParScp(5); pB = pBx0;
Nx = mat_getParScp(6); Np = mat_getParScp(7);
Dt = mat_getParScp(8); numTimeSteps = mat_getParScp(9);
cMovieFrameFreq = mat_getParScp(10);
if (length(mat_getParScp) >= 11)  % Backward compatible with old formats
    aveFreq = mat_getParScp(11);
end

% Calculate related parameters
Dx = (xf - x0) / Nx; Dp = (pB - pA) / Np;
numCellsX = Nx + 2; numCellsP = Np + 2;

clear mat_getParScp