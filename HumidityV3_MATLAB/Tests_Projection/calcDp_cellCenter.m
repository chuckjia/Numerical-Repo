function Dp_vec = calcDp_cellCenter(x0, pA, Dx, Nx, Np, cellCentersX)
%CALCDP_CELLCENTER Calculate Dp, i.e. the p step size, for all cells
%   OUTPUT:: DpVec: A vector of size Nx, whose ith entry is the Dp value at the center of the (i,j) cell

% % Approach 1
% Dp_vec = zeros(Nx);
% for i = 1:Nx
%     xLeft = x0 + (i - 1) * Dx;
%     xRight = xLeft + Dx;
%     cellLeftDp = (pB_fcn(xLeft) - pA) / Np;
%     cellRightDp = (pB_fcn(xRight) - pA) / Np;
%     Dp_vec(i) = 0.5 * (cellLeftDp + cellRightDp);
% end

% Approach 2
xVec = cellCentersX(:,1);
Dp_vec = (pB_fcn(xVec) - pA) ./ Nx;

end