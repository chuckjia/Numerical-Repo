 folder = "Storage/Res11/"
% folder = "Storage/Res7_Nx400_Dt0.125/"
% folder = "Results/"
%tVec = 50:50:5000;
 tVec = getVecFromFile_fcn(folder, "time");
norm_T = getVecFromFile_fcn(folder, "T_norm");
dotSize = 10;
figure
plot(tVec, norm_T, '.-', 'MarkerSize', dotSize); title("T")

norm_q = getVecFromFile_fcn(folder, "q_norm");
figure
plot(tVec, norm_q, '.-', 'MarkerSize', dotSize); title("q")

norm_u = getVecFromFile_fcn(folder, "u_norm");
figure
plot(tVec, norm_u, '.-', 'MarkerSize', dotSize); title("u")

norm_w = getVecFromFile_fcn(folder, "w_norm");
figure
plot(tVec, norm_w, '.-', 'MarkerSize', dotSize); title("w")
