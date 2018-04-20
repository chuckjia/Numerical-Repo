clear; clc

path = "/Users/chuckjia/Documents/Workspace/DataStorage/Humidity/";
folder = path + "Sim14/Norm/";
dotSize = 10;

graph_T = 1;
graph_q = 1;
graph_u = 1;
graph_w = 1;
saveAsPDF = true;

tVec = getVecFromFile_fcn(folder, "time");

if (graph_T)
    startLevel_T = 1.987e6;  % 2.055e6 for original simulation
    norm_T = getVecFromFile_fcn(folder, "T_norm");
    figure; plot(tVec, norm_T - startLevel_T, '.-', 'MarkerSize', dotSize)
    title("T"); xlabel("Time"); ylabel(strcat('Value - ', num2str(startLevel_T, '%1.4e')))
    if (saveAsPDF) print('Results/norm_T', '-dpdf', '-bestfit'); end
end

if (graph_q)
    norm_q = getVecFromFile_fcn(folder, "q_norm");
    figure; plot(tVec, norm_q, '.-', 'MarkerSize', dotSize);
    title("q"); xlabel("Time"); ylabel("Value")
    if (saveAsPDF) print('Results/norm_q', '-dpdf', '-bestfit'); end
end

if (graph_u)
    norm_u = getVecFromFile_fcn(folder, "u_norm");
    figure; plot(tVec, norm_u, '.-', 'MarkerSize', dotSize);
    title("u"); xlabel("Time"); ylabel("Value")
    if (saveAsPDF) print('Results/norm_u', '-dpdf', '-bestfit'); end
end

if (graph_w)
    norm_w = getVecFromFile_fcn(folder, "w_norm");
    figure; plot(tVec, norm_w, '.-', 'MarkerSize', dotSize);
    title("w"); xlabel("Time"); ylabel("Value")
    if (saveAsPDF) print('Results/norm_w', '-dpdf', '-bestfit'); end
end

