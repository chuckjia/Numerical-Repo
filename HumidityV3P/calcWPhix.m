function [w_, phix_] = calcWPhix(param, T_, u_, w_, phix_, cellCentersP, a1, a2mat, a3mat, a4mat, M_e12, M_e22)
%CALCW Summary of this function goes here
%   Detailed explanation goes here

Dx = param.Dx;
Np = param.Np;

gradxU = calcGradx(Dx, u_, a1, a2mat, a3mat, a4mat, M_e12, M_e22);
w_(:, 1) = 0;

for j = 0:Np-1
    jj = j + 1;
    p_factor = cellCentersP(:, jj+1) - cellCentersP(:, jj);
    w_(:, jj+1) = w_(:, jj) - [p_factor(1:end-1) .* gradxU(:, jj); 0];
end

% gradxT = calcGradx(Dx, T_, a1, a2mat, a3mat, a4mat, M_e12, M_e22);


end

