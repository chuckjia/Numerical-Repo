function y = pB_fcn(x)
%PB_FCN pB function as in the paper, describing the mountain surface

y = 1000 - 250 .* exp(-((x-37500) ./ 6000).^2);

end

