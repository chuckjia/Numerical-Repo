function intU_vec = calcIntU(u_mat, Dp_vec)
%CALCINTU Calculate the integral of u with respect to p, i.e. int_pA^pB u(x,p) dp.
%   INPUTS:: u_mat: Matrix of size Nx by Np. Its (i,j) entry is the numerical value of u at (i,j) cell center
%            Dp_mat: Matrix of size Nx by Np. Its (i,j) entry is the Dp value at the center of the (i,j) cell
%   OUTPUT:: A vector of size Nx, whose ith entry is the value int_pA^pB u(x_i,p) dp. Here x_i is the
%              x-coordinate of the center of the (i,j) cell, j = 1, ..., Np

Nx = length(u_mat);
intU_vec = zeros(Nx, 1);
for i = 1:Nx
    intU_vec(i) = sum(u_mat(i, :)) .* Dp_vec(i);
end

end