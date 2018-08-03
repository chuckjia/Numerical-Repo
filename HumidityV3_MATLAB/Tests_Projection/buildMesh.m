function [meshGridX, meshGridP] = buildMesh(x0, xf, pA, Nx, Np) 
%BUILDMESH Calculate the coordinate of the mesh grids
%   OUTPUT:: meshGridX: a matrix of size Nx+1 by Np+1, whose (i,j) entry is x_{i-1/2}
%            meshGridP: a matrix of size Nx+1 by Np+1, whose (i,j) entry is p_{i-1/2, j-1/2}

meshGridX = zeros(Nx + 1, Np + 1);  % x-coordinates of the mesh grids
meshGridP = zeros(Nx + 1, Np + 1);  % p-coordinates of the mesh grids
Dx = (xf - x0) / Nx;

for i = 1:Nx+1
    cellLeftX = x0 + (i - 1) * Dx;
    Dp = (pB_fcn(cellLeftX) - pA) / Np;
    for j = 1:Np+1
        meshGridX(i, j) = cellLeftX;
        meshGridP(i, j) = pA + (j - 1) * Dp;
    end
end

end

