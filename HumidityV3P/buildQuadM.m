function [M_e12, M_e22] = buildQuadM(meshGridP, cellCentersP)
%BUILDQUADMINV Calculates the first row of the inverse of the matrix M_{i,j+1/2} as in P108

Nx = size(cellCentersP, 1) - 2;
Np = size(cellCentersP, 2) - 2;

M_e12 = zeros(Nx+1, Np+1);
M_e22 = zeros(Nx+1, Np+1);

for ii = 1:Nx+1
    for jj = 1:Np+1
        M_e12(ii, jj) = meshGridP(ii+1, jj+1) - meshGridP(ii, jj+1);
        M_e22(ii, jj) = cellCentersP(ii, jj+1) - cellCentersP(ii, jj);
    end
end

end

