clearAllScp

x0 = 0;  xf = 75000; pA = 250;
meshN = 400;  Nx = meshN; Np = meshN;
Dx = (xf-x0) / Nx;
meshX = x0:Dx:xf;

outputFolder = '~/git/Numerical-Repo/HumidityV3/Output/';
cellCenterX = csvread(outputFolder + "CellCenters_X.csv");  cellCenterP = csvread(outputFolder + "CellCenters_P.csv");
cellCenterX = cellCenterX(2:end-1, 2);  cellCenterP = cellCenterP(2:end-1, 2:end-1);
DpMat = calc_DpMat(Nx, Np, cellCenterX, x0, pA, Dx);

a_vec = pBx_fcn(cellCenterX);
b_vec = pB_fcn(cellCenterX) - pA;
uTildeMat = uTilde_fcn(cellCenterX, cellCenterP);
int_u_vec = calc_int_u_vec(uTildeMat, DpMat);

c_vec = zeros(Nx, 1);
for i = 1:(Nx-1)
    c_vec(i) = (1/Dx) * (int_u_vec(i + 1) - int_u_vec(i));
end

A = zeros(Nx);

for i = 1:(Nx-1)
    A(i, i) = a_vec(i) - b_vec(i)/Dx;
    A(i, i+1) = b_vec(i)/Dx;
end
A(Nx, :) = 1;

lambdax_vec = A\c_vec;

c_res = csvread(outputFolder + "lambdax_diagnostics.csv")';
lambdax_vec - c_res

uMat_after_proj = apply_proj(uTildeMat, c_res);
new_int_u_vec = calc_int_u_vec(uMat_after_proj, DpMat);
max(abs(calc_der_x_int_u(new_int_u_vec, Dx)))


%%

function y = pB_fcn(x)

y = 1000 - 250 .* exp(-((x-37500) ./ 6000).^2);

end

function y = pBx_fcn(x)

y = (1/72000) .* (x - 37500) .* exp(-((x-37500) ./ 6000).^2);

end

function y = uTilde_fcn(x, p)

p0 = 1000;
xf = 75000;
y = 7.5 + 2 .* cos(p .* pi ./ p0) .* cos(2 .* pi .* x ./ xf);

end

function DpMat = calc_DpMat(Nx, Np, cellCenterX, x0, pA, Dx)
DpMat = zeros(Nx, Np);
for i = 1:Nx
    xCenter = cellCenterX(i);
    xLeft = x0 + (i - 1) * Dx;
    xRight = xLeft + Dx;
    cellLeftDp = (pB_fcn(xLeft) - pA) / Np;
    cellRightDp = (pB_fcn(xRight) - pA) / Np;
    r = (xCenter - xLeft) / Dx;
    
    for j = 1:Np
        pBottLeft = pA + (j - 1) * cellLeftDp;
        pTopLeft = pA + j * cellLeftDp;
        pBottRight = pA + (j - 1) * cellRightDp;
        pTopRight = pA + j * cellRightDp;
        
        pTopCenter = pTopLeft + r * (pTopRight - pTopLeft);
        pBottCenter = pBottLeft + r * (pBottRight - pBottLeft);
        
        DpMat(i, j) = pTopCenter - pBottCenter;
    end
end
end

function res = calc_der_x_int_u(int_u_ver, Dx)

res_len = length(int_u_ver) - 1;
res = zeros(res_len, 1);
for i = 1:res_len
    res(i) = (int_u_ver(i + 1) - int_u_ver(i)) / Dx;
end

end

function u = apply_proj(u, lambdax)

Nx = length(lambdax);
for i = 1:Nx
    u(i, :) = u(i, :) - lambdax(i);
end

end

function int_u_vec = calc_int_u_vec(uMat, DpMat)

Nx = length(uMat);
int_u_vec = zeros(Nx, 1);
for i = 1:Nx
    int_u_vec(i) = sum(uMat(i, :) .* DpMat(i, :));
end

end














