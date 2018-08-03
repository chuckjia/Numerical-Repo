function y = uTilde_fcn(x, p)
% UTILDE_FCN The initial condition for u

p0 = 1000;
xf = 75000;
y = 7.5 + 2 .* cos(p .* pi ./ p0) .* cos(2 .* pi .* x ./ xf);

end