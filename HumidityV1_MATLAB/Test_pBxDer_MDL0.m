clear; clc

syms x

pB_fcn = 1000. - 250. * exp(- (x - 37500.)^2 / 6000^2);
pB_xDer_fcn = diff(pB_fcn, x);

arr = [1, 3e4, 3.5e4, 4e4, 4.5e4, 5e4, 6e4, 7e4];
n = length(arr);

for i = 1:n
    x = arr(i);
    fprintf("%1.7e\n", double(subs(pB_xDer_fcn)));
end

