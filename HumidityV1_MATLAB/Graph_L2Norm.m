folder = "Storage/Res17/"
dotSize = 10;

tVec = getVecFromFile_fcn(folder, "time");

norm_T = getVecFromFile_fcn(folder, "T_norm");
figure
plot(tVec, norm_T - 2.055e6, '.-', 'MarkerSize', dotSize)
title("T"); xlabel("Time"); ylabel("Value - 2.055e6")

norm_q = getVecFromFile_fcn(folder, "q_norm");
figure
plot(tVec, norm_q, '.-', 'MarkerSize', dotSize); 
title("q"); xlabel("Time"); ylabel("Value")

norm_u = getVecFromFile_fcn(folder, "u_norm");
figure
plot(tVec, norm_u, '.-', 'MarkerSize', dotSize); 
title("u"); xlabel("Time"); ylabel("Value")

norm_w = getVecFromFile_fcn(folder, "w_norm");
figure
plot(tVec, norm_w, '.-', 'MarkerSize', dotSize); 
title("w"); xlabel("Time"); ylabel("Value")
