clear; clc; close all
projectPath = "~/Documents/Workspace/git/Numerical-Repo/HumidityV3/";

% Mesh grid: cell centers
param = readParam(projectPath + "Output/Param.csv");
x0 = param.x0;  xf = param.xf;  Nx = param.Nx;  Np = param.Np;  Dx = (xf - x0) / Nx;
centersX = csvread(projectPath + "Output/CellCenters_X.csv");
centersP = csvread(projectPath + "Output/CellCenters_P.csv");

solnName = 'q';
stepNo = -1;
plotName = "Numerical Solution";
plotLevel = Np - 60;

resultFolder = "MovieFrames/";
stepNo = genActualStepNo(stepNo, projectPath + resultFolder, solnName);

% Generate graph titles
titleLine1 = "Plot of " + plotName + " " + solnName;
plotTime = param.Dt * stepNo;
titleLine2 = "Time = " + num2str(plotTime) + "s";
viewAngle = [0, -90];  % Viewing angle
titleLine3 = "j = " + num2str(plotLevel) + ", Np = " + num2str(Np);

solnFileFullPath = genSolnFileFullPath(projectPath, resultFolder, solnName, stepNo);
[~, soln] = graphSoln(solnFileFullPath, centersX, centersP, titleLine1, titleLine2, viewAngle);

figure
x_vec = centersX(:, plotLevel);
soln_vec = soln(:, plotLevel);
plot(x_vec, soln_vec, 'LineWidth', 2);
max_soln = max(soln_vec);  min_soln = min(soln_vec);
unit_height = max_soln - min_soln;

hold on
pB = @(x) 1000 - 250 * exp(-((x - 37500) / 6000) .^ 2);
pB_vec = pB(x_vec);
max_pB = max(pB_vec);
magnify = unit_height;  shiftUp = max_soln - magnify + (max_soln - min_soln) * 0.5;
pB_vec = pB_vec ./ max_pB * magnify + shiftUp;
fig = plot(x_vec, pB_vec, '-.', 'Color', 0.75*[1,1,1], 'LineWidth', 1);
xlabel('x-axis');  ylabel('Value of Solution q');
legend('Solution q', 'Mountain (for reference)','Location', 'southwest');
title({titleLine1, titleLine2, titleLine3, ""})
hold off
filename = "Output/" + solnName + "_j_" + num2str(plotLevel);
saveas(fig, char(filename), 'pdf');

