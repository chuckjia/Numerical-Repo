function main()
%MAIN Build mesh grid and apply projection method on the initial u values

% ===== ===== ===== ===== ===== ===== 
% Mesh Paramters
% ===== ===== ===== ===== ===== ===== 

% Mountain geometry
x0 = 0;  xf = 75000;  pA = 250;

% Size of mesh used in the Finite Volume Method
Nx = 200;  Np = 200;

% Build mesh and calculate cell centers
Dx = (xf - x0) / Nx;  % x step size
[meshGridX, meshGridP] = buildMesh(x0, xf, pA, Nx, Np);  % meshGridX/meshGridP are the x-/p-coordinates of the mesh grid
[cellCentersX, cellCentersP] = calcCellCenters(meshGridX, meshGridP);  % cellCentersX/cellCentersP are the x-/p-coordinates of the cell centers

% Calculate Dp, i.e. the p step size, for all cells. The (i,j) entry of the output matrix is the Dp value at 
% the center of the cell
Dp_vec = calcDp_cellCenter(x0, pA, Dx, Nx, Np);

% ===== ===== ===== ===== ===== ===== 
% Projection Method
% ===== ===== ===== ===== ===== ===== 

% Now calculate cooefficients a_i, b_i, and c_i as in (3.34)
a_vec = pBx_fcn(cellCentersX, Dx);  % Coefficients a_i
b_vec = pB_fcn(cellCentersX) - pA;  % Coefficient b_i

uTilde_mat = uTilde_fcn(cellCentersX, cellCentersP);  % Initial condition of u
intU_vec = calcIntU(uTilde_mat, Dp_vec);  % Integral of u with respect to p, i.e. int_pA^pB u(x,p) dp
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
intU_afterProj = calcIntU(u_mat_afterProj, Dp_vec);  % int_pA^pB u dp after projection
result = calcDxIntU(intU_afterProj, Dx);  % (d/dx)int_pA^pB u dp at each x_i
[maxVal, maxInd] = max(abs(result));
fprintf("The max of |(d/dx)int_pA^pB u dp| is %1.4e, which occurs at the point x_%d\n", maxVal, maxInd);

end