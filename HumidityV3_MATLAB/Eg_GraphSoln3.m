
%% This script plots and updates graphs for solutions whenever new graph files are generated

clear; clc; close all

% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== 
% Settings
% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== 

% projectPath = "~/Documents/Workspace/Git/Numerical-Repo/HumidityV3/";
projectPath = "~/git/Numerical-Repo/HumidityV3/" % On Ubuntu
monitor = "left";

% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== 
% Parameter Setting and Validation: Do Not Modify
% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== 

currNumFiles = 0;
prevNumFiles = -1;
imageNo = 0;
itrNo = 0;
figPositionVec = [325, 500, 650, 500];  % [left bottom width height]
[accname, accno] = whichcomputer();
if accno == 1
    monitor = "right";
end
if monitor == "left"
    figPositionVec = [-1750, 500, 650, 500];  % [left bottom width height]
end
pauseTime1 = 1;
pauseTime2 = 1;

% ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
% Never-ending loop
% ===== ===== ===== ===== ===== ===== ===== ===== ===== =====

while true
    itrNo = itrNo + 1;
    if mod(itrNo, 2) == 0
        fprintf("... ")
    end
    prevNumFiles = currNumFiles;
    currNumFiles = countSolnFiles(projectPath);
    if (currNumFiles > prevNumFiles)
        close all
        figure('Position', figPositionVec)  
        imageNo = imageNo + 1;
        fprintf("Plotting graph no. %d.\n", imageNo);
        graphSolnRep("T", -1);
        figure('Position', figPositionVec + [figPositionVec(3)+25, 0, 0, 0])
        graphSolnRep("q", -1);
        % myCountdown(waittime)
        pause(pauseTime1)
    else
        pause(pauseTime2)
    end
    
end

%%

function myCountdown(n)

if nargin == 0 n = 30; end
fprintf("Countdown = ");
for j = 1:n
    fprintf("%5ds\n", n - j);
    pause(1)
    fprintf("\b\b\b\b\b\b\b");
end
fprintf("Updating now.\n");

end

%%

function n = countSolnFiles(projectPath)

n = length( dir(genFolderPathName(projectPath) + "MovieFrames") );

end

%%

function graphSolnRep(solnName, stepNo, projectPath)
% This file is an example of graphing surface and contour plots of numerical solutions.
%

if nargin < 3
    [~, accno] = whichcomputer();
    if accno == 1
        projectPath = "~/Documents/Workspace/Git/Numerical-Repo/HumidityV3/";
    else
        projectPath = "~/git/Numerical-Repo/HumidityV3/";  % Path on Ubuntu
    end
end
projectPath = genFolderPathName(projectPath);

% ===== ===== ===== =====
% Settings
% ===== ===== ===== =====

% solnName = "T";
% stepNo = -1;  % A value of -1 indicates graph the latest solution

resultFolder = "MovieFrames/";  % Commonly used: "MovieFrames/" or "Output/"
plotName = "Numerical Solution";  % "Solution", "Error", or "Exact Solution"

viewAngle = [0, -90];  % Viewing angle

% Setting some of the values: Do Not Change
stepNo = genActualStepNo(stepNo, projectPath + resultFolder, solnName, false);

% ===== ===== ===== =====
% Graph solution
% ===== ===== ===== =====

param = readParam(projectPath + "Output/Param.csv");
solnFileFullPath = genSolnFileFullPath(projectPath, resultFolder, solnName, stepNo);
centersX = csvread(projectPath + "Output/CellCenters_X.csv");
centersP = csvread(projectPath + "Output/CellCenters_P.csv");

% Generate graph titles
titleLine1 = "Plot of " + plotName + " " + solnName;
plotTime = param.Dt * stepNo;
titleLine2 = "Time = " + num2str(plotTime) + "s";
contourLevels = selectContourLevels(solnName);
contourProportion = 0.5;
graphContour(solnFileFullPath, centersX(2:end-1, 2:end-1), centersP(2:end-1, 2:end-1), contourLevels, titleLine1, titleLine2, ...
    viewAngle, contourProportion);

end