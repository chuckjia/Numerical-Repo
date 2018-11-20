
%% 

clear; clc; close all

solnName = 'T';
stepNo = -1;
% projectPath = "~/Documents/Workspace/git/Numerical-Repo/HumidityV3/";
projectPath = "/Users/chuckjia/Documents/Workspace/DataStorage/Humidity/HumidityV3/2018_11_6_1_14_14_h270_2700m/";


% iVec = [11:19];
iVec = [1:10, 20:10:50, 60:20:190, 195, 197, 199];
% iVec = [1:10, 20:10:50, 60:30:360, 370:10:390, 395];
% iVec = [1:20, 30:10:50, 60:50:760, 770:10:790, 795];
graphAtFixedP_wrapper(solnName, iVec, projectPath, stepNo);


%% 

function graphAtFixedP_wrapper(solnName, iVec, projectPath, stepNo) 

fprintf("Plotting %d levels in total.\n", length(iVec));
for j = 1:length(iVec)
    i = iVec(j);
    fprintf("Plotting at i = %d. Progress = %1.2f%%\n", i, j / length(iVec) * 100);
    graphAtFixedP_fcn(solnName, i, projectPath, stepNo);
end

end

%%

function graphAtFixedP_fcn(solnName, i, projectPath, stepNo)
% graphAtFixedP_fcn Generate plots of solutions at a fixed p level.
%    INPUT:: solnName    : The name of the solution. Currently, 2 options are supported: 'T' or 'q', with 'q'
%                              as the default
%            i           : The number of levels away from mountain surface. In other words, we are plotting  
%                              the level Np-i
%            projectPath : The path to the project directory
% 

close all
if nargin < 4
    stepNo = -1;
    if nargin < 3
        projectPath = "~/Documents/Workspace/git/Numerical-Repo/HumidityV3/";
        if nargin < 2
            i = 1;
            if nargin < 1
                solnName = 'q';
            end
        end
    end
end
projectPath = genFolderPathName(projectPath);

% Mesh grid: cell centers
param = readParam(projectPath + "Output/Param.csv");
x0 = param.x0;  xf = param.xf;  Nx = param.Nx;  Np = param.Np;  Dx = (xf - x0) / Nx;
centersX = csvread(projectPath + "Output/CellCenters_X.csv");
centersP = csvread(projectPath + "Output/CellCenters_P.csv");

%solnName = 'T';
% stepNo = -1;
plotName = "Numerical Solution";
plotLevel = Np - i;

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
max_pB = max(pB_vec);
pB_vec = pB_vec ./ max_pB * magnify + shiftUp;
max_pB = max(pB_vec);  yLimUnit = (max_pB - min_soln) * 0.2;
fig = plot(x_vec, pB_vec, '-', 'Color', 0.75*[1,1,1], 'LineWidth', 1);
ylim([min_soln - yLimUnit, max_pB + yLimUnit]);
xlabel('x-axis');  ylabel('Value of Solution q');
title({titleLine1, titleLine2, titleLine3, ""})
x_middle = (x0 + xf) * 0.5;
plot([x_middle, x_middle], [min_soln - yLimUnit, max_pB + yLimUnit], '--', 'Color', 0.9*[1,1,1], 'LineWidth', 0.5)
legend('Solution q', 'Mountain (for reference)','Location', 'southwest');
hold off
filename = "Output/" + solnName + "_j_" + num2str(plotLevel);
saveas(fig, char(filename), 'pdf');


end


