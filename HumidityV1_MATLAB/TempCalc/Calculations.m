

x = 37500;
pB = pB_fcn(x);

R_CONST = 900.;  % 287.;
T0 = 300.; 
DT = 50.;
p0 = 1000.;

zB = (- R_CONST * (T0 - DT) * (log(pB) - log(p0)) - R_CONST * DT / p0 * pB + R_CONST * DT) / 9.8;

vpa(zB, 6)

x = 0:1000:75000;
y = pB_fcn(x);
plot(x, y)





















