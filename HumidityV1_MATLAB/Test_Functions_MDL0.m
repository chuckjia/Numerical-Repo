clear; clc
syms x
heightFactor = 250;

pB_fcn = 1000 - heightFactor * exp(-(x - 37500)^2 / 6000^2);
pB_xDer_fcn = diff(pB_fcn, x);

x = 37400;
fprintf("At x = %1.1f, pB = %1.10f, pB_x = %1.10f\n", ...
    x, double(subs(pB_fcn)), double(subs(pB_xDer_fcn)));
