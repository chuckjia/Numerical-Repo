clear; clc

folder = "Results/";
resFromC = getVecFromFile_fcn(folder, "test_exactSoln_MDL1");
clear folder

syms x p t pA xf R p0 T0 DeltaT

pB_fcn = 1000 - 200 * exp(- (x - 25000)^2 / 9e6);
xi_fcn = (p - pA)^3 * (p - pB_fcn)^3 * 1e-12 * (cos(2 * pi * t) + 20) * x^3 * (x - xf)^3 / xf^6;
phi_fcn = ( ...
    (p - pB_fcn)^3 / 450^3 - R / 9.8 * ( (T0 - DeltaT) * log(p) + DeltaT / p0 * p) ...
    ) * cos(2 * pi * t) * x * (x - xf)^2 / xf^3;

T_exact = -p / R * diff(phi_fcn, p);
T_xDer_exact = diff(T_exact, x);
T_pDer_exact = diff(T_exact, p);
T_tDer_exact = diff(T_exact, t);
u_exact = -diff(xi_fcn, p);
u_xDer_exact = diff(u_exact, x);
u_pDer_exact = diff(u_exact, p);
u_tDer_exact = diff(u_exact, t);
w_exact = diff(xi_fcn, x);
w_pDer_exact = diff(w_exact, p);

pA = 200; xf = 50000; R = 287; p0 = 1000; T0 = 300; DeltaT = 50;
x = resFromC(1); p = resFromC(2); t = resFromC(3);

res = double([
    subs(T_exact);
    subs(T_xDer_exact); subs(T_pDer_exact); subs(T_tDer_exact);
    subs(u_exact); 
    subs(u_xDer_exact); subs(u_pDer_exact); subs(u_tDer_exact);
    subs(w_exact); subs(w_pDer_exact)]);

fprintf("Evaluating at the point x = %1.2f, p = %1.2f, t = %1.2fs\n", ... 
    double(x), double(p), double(t));
err = abs(resFromC(4:13) - res);
maxErr = max(err);
fprintf("The largest error is %1.2e\n", maxErr); 
vpa(err.', 9)






