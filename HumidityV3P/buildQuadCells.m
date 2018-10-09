function [a2mat, a3mat, a4mat] = buildQuadCells(meshGridX, meshGridP, cellCentersX, cellCentersP, a1)
%BUILDQUADCELLS This function calculates and stores the quadrilateral cell coefficients

Nx = size(cellCentersX, 1) - 2;
Np = size(cellCentersX, 2) - 2;

a2mat = zeros(Nx+1, Nx+1);
a3mat = zeros(Nx+1, Nx+1);
a4mat = zeros(Nx+1, Nx+1);


for ii = 1:Nx+1
    for jj = 1:Np+1        
        x1 = cellCentersX(ii, jj);
        x2 = cellCentersX(ii+1, jj);
        x3 = cellCentersX(ii, jj+1);
        x4 = cellCentersX(ii+1, jj+1);
        
        p1 = cellCentersP(ii, jj);
        p2 = cellCentersP(ii+1, jj);
        p3 = cellCentersP(ii, jj+1);
        p4 = cellCentersP(ii+1, jj+1);
        
        x_rhs = meshGridX(ii+1, jj+1);
        p_rhs = meshGridP(ii+1, jj+1);
        
        A = [x2, x3, x4;
            p2, p3, p4;
            1,  1,  1];
        b = [x_rhs - a1 * x1; p_rhs - a1 * p1; 1 - a1];
        a234 = A\b;
        
        a2mat(ii, jj) = a234(1);
        a3mat(ii, jj) = a234(2);
        a4mat(ii, jj) = a234(3);
    end
end

