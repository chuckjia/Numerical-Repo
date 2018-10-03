function y = pBx_fcn(cellCentersX, Dx)
%PBX_FCN The x-derivative of the pB function

% len = length(cellCentersX);
% y = zeros(len, 1);
% for i = 1:(len-1)
%     y(i) = (pB_fcn(cellCentersX(i + 1)) - pB_fcn(cellCentersX(i))) / Dx;
% end

term = (cellCentersX - 37500) ./ 6000;
y = exp(-term.^2) .* term ./ 12;

end