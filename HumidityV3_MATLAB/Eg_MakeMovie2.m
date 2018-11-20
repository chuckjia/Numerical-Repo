% This file is an example of making a movie of the numerical errors.
%

clearAllScp

% ===== ===== ===== ===== 
% Settings
% ===== ===== ===== ===== 

solnName = "T";
steps = 1:100;  % Vector of all the step numbers to be included in movie
zTopVal = 7e-5;
figureSize = [10, 10, 900, 700];  % x-pos, y-pos, width, height
% Path and solution file names
projectPath = "~/Documents/Workspace/Git/Numerical-Repo/HumidityV3/";  % Path to the outermost folder
outputFilename = "Output/Error_" + solnName + "_slow.avi";

% ===== ===== ===== ===== 
% Generate Movies
% ===== ===== ===== ===== 

% Mesh and paramters
param = readParam(projectPath + "Output/Param.csv"); 
centersX = csvread(projectPath + "Output/CellCenters_X.csv"); centersX = centersX(2:end-1, 2:end-1);
centersP = csvread(projectPath + "Output/CellCenters_P.csv"); centersP = centersP(2:end-1, 2:end-1);

titleLine1 = "Numerical Errors on " + solnName;

figure('pos', figureSize)
frameNo = 0;
for stepNo = steps
    frameNo = frameNo + 1;
    time = stepNo * param.Dt;
    titleLine2 = "time = " + num2str(time) + "s";
    filename = projectPath + "MovieFrames/err_" + solnName + "_" + int2str(stepNo) + ".csv";
    
    viewAngle = false;
    graphSoln(filename, centersX, centersP, titleLine1, titleLine2, viewAngle);
    zlim([-zTopVal, zTopVal])

    F(frameNo) = getframe(gcf);
    if (frameNo == 1) % Allocate memory space at the beginning
        F = repmat(F, 1, length(steps));
    end
end

fprintf("Parameters used: mesh size = %dx%d\n", param.Nx, param.Np);
slowMotion(F, 5, char(outputFilename));  




