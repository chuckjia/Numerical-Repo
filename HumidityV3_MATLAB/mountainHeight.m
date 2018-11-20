function height = mountainHeight(h, x0, xf)
%MOUNTAINHEIGHT Summary of this function goes here
%   Detailed explanation goes here

x = 0.5 * (x0 + xf);
term = x - 37500.;
height = 1000 - h * exp(-term * term / 3.6e7);
height = p2z(height);

end

function z = p2z(p)

R = 287;
T0 = 300;
DeltaT = 50; 
p0 = 1000;

z = 1 / 9.8 * (R * (T0 - DeltaT) * log(p0 / p) - R * DeltaT / p0 * p + R * DeltaT);

end