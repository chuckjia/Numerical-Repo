function gradxMat = calcGradx(Dx, u, a1, a2mat, a3mat, a4mat, M_e12, M_e22)
%CALCGRADX Calculate the (grad_x u_h) and (grad_x T_h) using the quadrilateral cell interpolations
%
%   INPUT::
%       Dx    : A scalar representing the x step size
%       u     : A matrix of size Nx+2 by Np+2, representing the values of a solution on the cell centers
%       a1    : Scalar, the 1st quadrilateral cell coefficients
%       a2mat : Matrix, the 2nd quadrilateral cell coefficients
%               a3mat, a4mat are similarly defined, as the 3rd and 4th coefficients matrices
%       M_e12 : The (1,2)-entry in the matrix M_{i,j+1/2} as in P108
%               M_e22 is similarly defined, storing the (2,2)-entry of the same matrices
%

% u_quad is a matrix of size Nx+1 by Np+1. Its (i,j)-entry stores u_{i+1/2, j+1/2}
u_quad = ...
    a1     .*  u(1:end-1, 1:end-1)  + ...
    a2mat  .*  u(2:end,   1:end-1)  + ...
    a3mat  .*  u(1:end-1, 2:end)    + ...
    a4mat  .*  u(2:end,   2:end);

% e1mat and e2mat stores the first and the second element of the vector in (3.15)
e1mat = u_quad(2:end, :) - u_quad(1:end-1, :);
e1mat = [zeros(1, size(e1mat, 2)); e1mat];
e2mat = u(1:end-1, 2:end) - u(1:end-1, 1:end-1);

gradxMat = (e1mat - M_e12 .* e2mat ./ M_e22) ./ Dx;

end

