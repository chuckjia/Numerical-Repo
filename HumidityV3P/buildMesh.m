function [meshGridX, meshGridP, DpVec] = buildMesh(param) 
%BUILDMESH Calculate the coordinate of the mesh grids
%
%   OUTPUT:: 
%       meshGridX :  a matrix of size Nx+3 by Np+3, whose (i+1,j+1) entry is x_{i-1/2}
%       meshGridP :  a matrix of size Nx+3 by Np+3, whose (i+1,j+1) entry is p_{i-1/2, j-1/2}
%       DpVec     :  a vector of size Nx+3, whose (i+1) entry is Dp_{i-1/2}   
%

x0 = param.x0;
xf = param.xf;
pA = param.pA;
Nx = param.Nx;
Np = param.Np;

numGridPtsX = Nx + 3;
numGridPtsP = Np + 3;

meshGridX = zeros(numGridPtsX, numGridPtsP);  % x-coordinates of the mesh grids
meshGridP = zeros(numGridPtsX, numGridPtsP);  % p-coordinates of the mesh grids
DpVec = zeros(numGridPtsX);
Dx = (xf - x0) / Nx;

for ii = 1:numGridPtsX
    i = ii - 1;
    cellLeftX = x0 + (i - 1) * Dx;
    Dp = (pB_fcn(cellLeftX) - pA) / Np;
    DpVec(ii) = Dp;
    for jj = 1:numGridPtsP
        j = jj - 1;
        meshGridX(ii, jj) = cellLeftX;
        meshGridP(ii, jj) = pA + (j - 1) * Dp;
    end
end

end

