function y = pBx_fcn(cellCentersX, Dx)
%PBX_FCN The x-derivative of the pB function

len = length(cellCentersX);
y = zeros(len + 1);
for i = 1:len
    y(i) = (pB_fcn(cellCentersX(i + 1)) - pB_fcn(cellCentersX(i))) / Dx;
end

end