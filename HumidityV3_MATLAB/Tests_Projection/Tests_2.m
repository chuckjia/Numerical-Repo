clearAllScp

% ===== ===== ===== ===== ===== ===== 
% Mesh Paramters
% ===== ===== ===== ===== ===== ===== 

% Mountain geometry
x0 = 0;
xf = 75000;
pA = 250;

% Size of mesh used in the Finite Volume Method
Nx = 200;
Np = 200;

% Build mesh and calculate cell centers
Dx = (xf - x0) / Nx;  % x step size
[meshGridX, meshGridP] = buildMesh(x0, xf, pA, Nx, Np);
[cellCentersX, cellCentersP] = calcCellCenters(meshGridX, meshGridP);

% Calculate Dp, i.e. the p step size, for all cells. The (i,j) entry of the output matrix is the Dp value at 
% the center of the cell
Dp_mat = calcDp_cellCenter(x0, pA, Dx, Nx, Np, cellCentersX);

% ===== ===== ===== ===== ===== ===== 
% Projection Method
% ===== ===== ===== ===== ===== ===== 

% Now calculate cooefficients a_i, b_i, and c_i as in (3.34)
a_vec = pBx_fcn(cellCentersX);  % Coefficients a_i
b_vec = pB_fcn(cellCentersX) - pA;  % Coefficient b_i

uTilde_mat = uTilde_fcn(cellCentersX, cellCentersP);  % Initial condition of u
intU_vec = calcIntU(uTilde_mat, Dp_mat);  % Integral of u with respect to p, i.e. int_pA^pB u(x,p) dp
c_vec = zeros(Nx, 1);  % Coefficient c_i
for i = 1:(Nx-1)
    c_vec(i) = (intU_vec(i + 1) - intU_vec(i)) / Dx;
end

% A represents the coefficient matrix of the linear system in (3.33), i.e. the LHS of (3.33)
A = zeros(Nx, Nx);
for i = 1:(Nx-1)
    A(i, i) = a_vec(i) - b_vec(i) / Dx;
    A(i, i+1) = b_vec(i) / Dx;
end
A(Nx, :) = 1;  % From (3.35)

% Solve linear system for lambda_x
lambdax_vec = A \ c_vec;

% Apply the projection method
u_mat_afterProj = applyProj(uTilde_mat, lambdax_vec);  % Projected u
intU_afterProj = calcIntU(u_mat_afterProj, Dp_mat);
calcDxIntU(intU_afterProj, Dx)


%%

function y = pB_fcn(x)
%PB_FCN pB function as in the paper, describing the mountain surface
y = 1000 - 250 .* exp(-((x-37500) ./ 6000).^2);

end

function y = pBx_fcn(x)
%PBX_FCN The x-derivative of the pB function
y = (1/72000) .* (x - 37500) .* exp(-((x-37500) ./ 6000).^2);

end

function y = uTilde_fcn(x, p)
% UTILDE_FCN The initial condition for u

p0 = 1000;
xf = 75000;
y = 7.5 + 2 .* cos(p .* pi ./ p0) .* cos(2 .* pi .* x ./ xf);

end

function [meshGridX, meshGridP] = buildMesh(x0, xf, pA, Nx, Np) 
% BUILDMESH Calculate the coordinate of the mesh grids

meshGridX = zeros(Nx + 1, Np + 1);
meshGridP = zeros(Nx + 1, Np + 1);
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

function [cellCentersX, cellCentersP] = calcCellCenters(meshGridX, meshGridP)
% CALCCELLCENTERS Calculate the coordinates of the cell centers

Nx = size(meshGridX, 1) - 1;
Np = size(meshGridX, 2) - 1;
cellCentersX = zeros(Nx, Np);
cellCentersP = zeros(Nx, Np);

for i = 1:Nx
    for j = 1:Np
        x_left = meshGridX(i, j);
        x_right = meshGridX(i + 1, j);
        p_bottLeft = meshGridP(i, j);
        p_bottRight = meshGridP(i + 1, j);
        p_topLeft = meshGridP(i, j + 1);
        p_topRight = meshGridP(i + 1, j + 1);
        
        x1 = (2 * x_left + x_right) / 3;
        p1 = (p_bottLeft + p_bottRight + p_topLeft) / 3;
        x2 = (x_left + 2 * x_right) / 3;
        p2 = (p_bottLeft + p_bottRight + p_topRight) / 3;        
        x3 = (x_left + 2 * x_right) / 3;
        p3 = (p_bottRight + p_topLeft + p_topRight) / 3;
        x4 = (2 * x_left + x_right) / 3;
        p4 = (p_bottLeft + p_topLeft + p_topRight) / 3;
        
        [cellCentersX(i, j), cellCentersP(i, j)] = calcLineIntersection(x1, p1, x2, p2, x3, p3, x4, p4);
    end
end

end

function [x, p] = calcLineIntersection(x1, p1, x2, p2, x3, p3, x4, p4) 
%CALCLINEINTERSECTION Calculate coordinate of the intersection point of two lines. Each line is represented by 
%   2 points on that line. We assume that the points (x1, p1) and (x3, p3) are on line 1. The other 2 points
%   are on line 2

k1 = (p3 - p1) / (x3 - x1);
k2 = (p4 - p2) / (x4 - x2);
x = (k1 * x1 - p1 - k2 * x2 + p2) / (k1 - k2);
p = k1 * (x - x1) + p1;

end

function DpMat = calcDp_cellCenter(x0, pA, Dx, Nx, Np, cellCentersX)
%CALCDP_CELLCENTER Calculate Dp, i.e. the p step size, for all cells
%   OUTPUT:: DpMat: A matrix of size Nx by Np, whose (i,j) entry is the Dp value at the center of the (i,j) cell 

DpMat = zeros(Nx, Np);

% Dp calculated as the length of the vertical line section passing the cell centers
for i = 1:Nx
    xCenter = cellCentersX(i);
    xLeft = x0 + (i - 1) * Dx;
    xRight = xLeft + Dx;
    cellLeftDp = (pB_fcn(xLeft) - pA) / Np;
    cellRightDp = (pB_fcn(xRight) - pA) / Np;
    r = (xCenter - xLeft) / Dx;
    
    for j = 1:Np
        pBottLeft = pA + (j - 1) * cellLeftDp;
        pTopLeft = pA + j * cellLeftDp;
        pBottRight = pA + (j - 1) * cellRightDp;
        pTopRight = pA + j * cellRightDp;
        
        pTopCenter = pTopLeft + r * (pTopRight - pTopLeft);
        pBottCenter = pBottLeft + r * (pBottRight - pBottLeft);
        
        DpMat(i, j) = pTopCenter - pBottCenter;
    end
end

% % Dp calculated as the average of the Dp values on the cell left and right sides
% for i = 1:Nx
%     xLeft = x0 + (i - 1) * Dx;
%     xRight = xLeft + Dx;
%     cellLeftDp = (pB_fcn(xLeft) - pA) / Np;
%     cellRightDp = (pB_fcn(xRight) - pA) / Np;
%     
%     DpMat(i, :) = 0.5 * (cellLeftDp + cellRightDp);
% end

end

function intU_vec = calcIntU(u_mat, Dp_mat)
%CALCINTU Calculate the integral of u with respect to p, i.e. int_pA^pB u(x,p) dp.
%   OUTPUT:: A vector of size Nx, whose ith entry is the value int_pA^pB u(x_i,p) dp. Here x_i is the
%              x-coordinate of the center of the (i,j) cell, j = 1, ..., Np

Nx = length(u_mat);
intU_vec = zeros(Nx, 1);
for i = 1:Nx
    intU_vec(i) = sum(u_mat(i, :) .* Dp_mat(i, :));
end

end

function res_vec = calcDxIntU(intU_vec, Dx)
%CALCDXINTU Calculate the x-derivative of the int_pA^pB u dp

resLen = length(intU_vec) - 1;
res_vec = zeros(resLen, 1);
for i = 1:resLen
    res_vec(i) = (intU_vec(i + 1) - intU_vec(i)) / Dx;
end

end

function u_mat = applyProj(u_mat, lambdax_vec)
% APPLYPROJ Apply the projection method on u

Nx = length(lambdax_vec);
for i = 1:Nx
    u_mat(i, :) = u_mat(i, :) - lambdax_vec(i);
end

end














