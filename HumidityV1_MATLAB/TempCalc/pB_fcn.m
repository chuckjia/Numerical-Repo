function y = pB_fcn( x )
%PB_FCN Summary of this function goes here
%   Detailed explanation goes here

    y = 1000. - 279. * exp(-((x - 37500.) / 6000.) .^ 2);

end

