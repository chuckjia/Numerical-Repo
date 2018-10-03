function [u_afterProj, u_beforeProj, derAfterProj, derBeforeProj, ...
    a_vec, b_vec, c_vec, intU_vec, lambdax_vec] = main(Nx)
%MAIN Build mesh grid and apply projection method on the initial u values
%    INPUT::
%             Nx: The number of spacial steps in the x-direction. The p-direction is assumed to have the same
%                 number of steps Np = Nx. The default value is 100, if no number is passed in.
%
%    OUTPUT:: 
%             u_afterProj: A matrix representing the values of u on the spacial mesh AFTER projection.
%             u_beforeProj: A matrix representing the values of u on the spacial mesh BEFORE projection.
%

% ===== ===== ===== ===== ===== ===== 
% Mesh Paramters
% ===== ===== ===== ===== ===== ===== 

% Size of mesh used in the Finite Volume Method
if nargin < 1  Nx = 100;  end  % Default value for numDivisions
Np = Nx;

% Mountain geometry
x0 = 0;  xf = 75000;  pA = 250;

% Build mesh and calculate cell centers
Dx = (xf - x0) / Nx;  % x step size
[meshGridX, meshGridP] = buildMesh(x0, xf, pA, Nx, Np);  % meshGridX/meshGridP are the x-/p-coordinates of the mesh grid
[cellCentersX, cellCentersP] = calcCellCenters(meshGridX, meshGridP);  % cellCentersX/cellCentersP are the x-/p-coordinates of the cell centers

% Calculate Dp, i.e. the p step size, for all cells. The (i,j) entry of the output matrix is the Dp value at 
% the center of the cell
Dp_vec = calcDp_cellCenter(x0, pA, Dx, Nx, Np, cellCentersX);

% ===== ===== ===== ===== ===== ===== 
% Projection Method
% ===== ===== ===== ===== ===== ===== 

% Now calculate cooefficients a_i, b_i, and c_i as in (3.34)
a_vec = pBx_fcn(cellCentersX(1:Nx-1, 1), Dx);  % Coefficients a_i
b_vec = pB_fcn(cellCentersX(1:Nx-1, 1)) - pA;  % Coefficient b_i

u_beforeProj = uTilde_fcn(cellCentersX, cellCentersP);  % Initial condition of u
intU_vec = calcIntU(u_beforeProj, Dp_vec);  % Integral of u with respect to p, i.e. int_pA^pB u(x,p) dp
% result = calcDxIntU(intU_vec, Dx);
% plot(cellCentersX(:,1), [result; 0]);

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
u_afterProj = applyProj(u_beforeProj, lambdax_vec);  % Projected u
fprintf("Applied the projection method on the initial u with mesh size %dx%d.\n", Nx, Np);

% Evaluation: before projection
intU_beforeProj = calcIntU(u_beforeProj, Dp_vec);  % int_pA^pB u dp after projection
derBeforeProj = calcDxIntU(intU_beforeProj, Dx);  % (d/dx)int_pA^pB u dp at each x_i
[maxVal, maxInd] = max(abs(derBeforeProj));
fprintf("Before projection: The max of |(d/dx)int_pA^pB u dp| is %1.4e, which occurs at the point x_%d=%1.4e.\n", ... 
    maxVal, maxInd, cellCentersX(maxInd, 1));

% Evaluation: after projection
intU_afterProj = calcIntU(u_afterProj, Dp_vec);  % int_pA^pB u dp after projection
derAfterProj = calcDxIntU(intU_afterProj, Dx);  % (d/dx)int_pA^pB u dp at each x_i
[maxVal, maxInd] = max(abs(derAfterProj));

fprintf("After projection: The max of |(d/dx)int_pA^pB u dp| is %1.4e, which occurs at the point x_%d=%1.4e.\n", ...
    maxVal, maxInd, cellCentersX(maxInd, 1));

end