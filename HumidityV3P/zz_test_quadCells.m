%% Tests on the quadrilateral cell methods
% This file contains all the tests designed for the quadrilateral cell interpolation

%% Testings 1

main

cOutputFolderPath = "~/Documents/Workspace/Git/Numerical-Repo/HumidityV3/Output/";
a2mat_c = csvread(fullfile(cOutputFolderPath, "a2_quad.csv"));
a3mat_c = csvread(fullfile(cOutputFolderPath, "a3_quad.csv"));
a4mat_c = csvread(fullfile(cOutputFolderPath, "a4_quad.csv"));

[maxVal2, maxInd2] = max(reshape(a2mat(2:end-1, 2:end-1) - a2mat_c(2:end-1, 2:end-1), [], 1));
[maxVal3, maxInd3] = max(reshape(a3mat(2:end-1, 2:end-1) - a3mat_c(2:end-1, 2:end-1), [], 1));
[maxVal4, maxInd4] = max(reshape(a4mat(2:end-1, 2:end-1) - a4mat_c(2:end-1, 2:end-1), [], 1));


%%

main
cOutputFolderPath = "~/Documents/Workspace/Git/Numerical-Repo/HumidityV3/Output/";

w_ = calcWPhix(param, T_, u_, w_, phix_, cellCentersP, a1, a2mat, a3mat, a4mat, M_e12, M_e22);
w_c = csvread(fullfile(cOutputFolderPath, "w_test.csv"));
max(reshape(w_(2:end-1, 2:end-1) - w_c(2:end-1, 2:end-1), [], 1))

surf(cellCentersX(2:end-1, 2:end-1), cellCentersP(2:end-1, 2:end-1), w_(2:end-1, 2:end-1))
figure
surf(cellCentersX(2:end-1, 2:end-1), cellCentersP(2:end-1, 2:end-1), w_c(2:end-1, 2:end-1))



%%





















