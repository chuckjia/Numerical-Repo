function param = getPar( file_path )
%GETPAR_SCP Summary of this function goes here
%   Detailed explanation goes here

mat = getVecFromFile(file_path);
param = ParSet;

param.x0 = mat(1);
param.xf = mat(2);
param.pA = mat(3);
param.pBx0 = mat(4);
param.pBxf = mat(5);
param.pB = param.pBx0;
param.Nx = mat(6);
param.Np = mat(7);
param.Dt = mat(8);
param.numTimeSteps = mat(9);
param.cMovieFrameFreq = mat(10);
param.aveFreq = mat(11);

param.Dx = (param.xf - param.x0) / param.Nx;
param.Dp = (param.pB - param.pA) / param.Np;
param.numCellsX = param.Nx + 2;
param.numCellsP = param.Np + 2;

end

