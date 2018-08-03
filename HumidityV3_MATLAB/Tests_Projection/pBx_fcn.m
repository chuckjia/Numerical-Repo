function y = pBx_fcn(x)
%PBX_FCN The x-derivative of the pB function

y = (1/72000) .* (x - 37500) .* exp(-((x-37500) ./ 6000).^2);

end