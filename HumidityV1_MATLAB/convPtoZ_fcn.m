function z = convPtoZ_fcn( p )
%CONVPTOZ_FCN Summary of this function goes here
%   Detailed explanation goes here
    R = 287; T0 = 300; DT = 50; p0 = 1000;
    z = (- R * (T0 - DT) * (log(p) - log(p0)) - R * DT / p0 * p + R * DT) / 9.8;
end

