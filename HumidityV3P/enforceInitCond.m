function [T_, q_, u_] = enforceInitCond(param, cellCentersX, cellCentersP, modelNo)
%ENFORCEINITCOND Summary of this function goes here
%   Detailed explanation goes here

if nargin < 4
    modelNo = 0;
end
if modelNo ~= 0
    modelNo = 0;
end

if modelNo == 0
    [T_, q_, u_] = initTqu_MDL0(cellCentersX, cellCentersP, param.xf);
end

end

function val = qs(T, p)
    val = 0.622 * 6.112 .* exp(17.67 * (T - 273.15) ./ (T - 29.65)) ./ p;
end

function [T_, q_, u_] = initTqu_MDL0(x, p, xf)

T0 = 300;
p0 = 1000;
DeltaT = 50;
n = 1;

T_ = T0 - (1 - p / p0) * DeltaT;
q_ = qs(T_, p) - 0.0052;
u_ = 7.5 + 2 * cos(pi / p0 * p) .* cos(2 * n * pi * x / xf);

end










