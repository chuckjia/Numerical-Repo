function param = setModelParam()
%SETMODELPARAM Assign values to model parameters

x0 = 0;
xf = 75000;
pA = 250;

% Mesh size
Nx = 200;
Np = 200;
Dt = 0.5;
Nt = 100;

% Other parameters
AveRate = 20;
MovieFrameRate = -1;
NumMsg = 100;

% Set model parameters
param = ModelParam(x0, xf, pA, Nx, Np, Dt, Nt, AveRate, MovieFrameRate, NumMsg);
% param.showParam();

end

