clear; clc
syms x

pB_fcn = 1000 - 250 * exp(-(x - 37500)^2 / 6000^2);
pB_xDer_fcn = diff(pB_fcn, x);

x = 35000;
vpa(double(subs(pB_xDer_fcn)), 10)