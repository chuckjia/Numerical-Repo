%% Tests on the quadrilateral cell methods
% This file contains all the tests designed for the quadrilateral cell interpolation

%% Testings 1: a1, ..., a4

main

cOutputFolderPath = "~/Documents/Workspace/Git/Numerical-Repo/HumidityV3/Output/";
a2mat_c = csvread(fullfile(cOutputFolderPath, "a2_quad.csv"));
a3mat_c = csvread(fullfile(cOutputFolderPath, "a3_quad.csv"));
a4mat_c = csvread(fullfile(cOutputFolderPath, "a4_quad.csv"));

[maxVal2, maxInd2] = max(reshape(a2mat(2:end-1, 2:end-1) - a2mat_c(2:end-1, 2:end-1), [], 1));
[maxVal3, maxInd3] = max(reshape(a3mat(2:end-1, 2:end-1) - a3mat_c(2:end-1, 2:end-1), [], 1));
[maxVal4, maxInd4] = max(reshape(a4mat(2:end-1, 2:end-1) - a4mat_c(2:end-1, 2:end-1), [], 1));

showMaxVal(maxVal2, "a2");
showMaxVal(maxVal3, "a3");
showMaxVal(maxVal4, "a4");

% Passed tests with n = 200, error ~ 1e-13

%% Testings 2: M matrix

main
e12_MInv_c = readcfile("MInv_e12_quad");
e12_MInv_m = -M_e12 ./ M_e22 / Dx;
maxVal = max(max(e12_MInv_c - e12_MInv_m));
showMaxVal(maxVal);

% Passed tests with n = 200, error ~ 1e-15

%% Testings 3: gradxU

main
u_quad_m = ...
    a1     .*  u_(1:end-1, 1:end-1)  + ...
    a2mat  .*  u_(2:end,   1:end-1)  + ...
    a3mat  .*  u_(1:end-1, 2:end)    + ...
    a4mat  .*  u_(2:end,   2:end);
u_quad_c = readcfile("u_quad");
u_quad_m(:, 2:end) - u_quad_c(:, 2:end)


%%
main
gradxU_c = readcfile("gradhUx");
gradxU_mat = calcGradx(Dx, u_, a1, a2mat, a3mat, a4mat, M_e12, M_e22);
% maxVal = max(max(gradxU_c(2:end, 2:end) - gradxU_mat(2:end, 2:end)));
% showMaxVal(maxVal, "gradxU");






%% Testings: w

main
cOutputFolderPath = "~/Documents/Workspace/Git/Numerical-Repo/HumidityV3/Output/";

w_ = calcWPhix(param, T_, u_, w_, phix_, cellCentersP, a1, a2mat, a3mat, a4mat, M_e12, M_e22);
w_c = csvread(fullfile(cOutputFolderPath, "w_test.csv"));
max(reshape(w_(2:end-1, 2:end-1) - w_c(2:end-1, 2:end-1), [], 1))

surf(cellCentersX(2:end-1, 2:end-1), cellCentersP(2:end-1, 2:end-1), w_(2:end-1, 2:end-1))
figure
surf(cellCentersX(2:end-1, 2:end-1), cellCentersP(2:end-1, 2:end-1), w_c(2:end-1, 2:end-1))



%%

function mat = readcfile(filename)

folderPath = "~/Documents/Workspace/Git/Numerical-Repo/HumidityV3/Output/";
filename = string(filename);
if ~endsWith(filename, '.csv')
   filename = filename + ".csv"; 
end
mat = csvread(fullfile(folderPath, filename));

end

function showMaxVal(maxVal, name) 
    if nargin < 2
        fprintf("The max value is %1.4e.\n", maxVal);
    else
        fprintf("The max value for %s is %1.4e.\n", name, maxVal);
    end
    
end

















