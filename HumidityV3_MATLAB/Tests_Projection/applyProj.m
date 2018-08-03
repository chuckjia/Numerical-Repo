function u_mat = applyProj(u_mat, lambdax_vec)
% APPLYPROJ Apply the projection method on u
%   INPUTS:: u_mat: Matrix of size Nx by Np. Its (i,j) entry is the numerical value of u at (i,j) cell center
%            lambdax_vec: Vector of length Nx. Its ith component is the alpha_i as in (3.32)
%   OUTPUT:: u_mat: Matrix of size Nx by Np, storing the values of projected u values

Nx = length(lambdax_vec);
for i = 1:Nx
    u_mat(i, :) = u_mat(i, :) - lambdax_vec(i);
end

end