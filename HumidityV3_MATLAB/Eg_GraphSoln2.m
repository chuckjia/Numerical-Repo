clear; clc

% Settings
solnName = "q";
stepNo = 10000;
saveToPDF = false;
viewAngle = [0, -90];

% Path and solution file names
path = "/Users/chuckjia/Documents/Workspace/Git/Numerical-Repo/HumidityV3/";  % Path to the outmost folder
solnFilename = genSolnFilename(path + "MovieFrames/", solnName, stepNo);

% Mesh and paramters
param = readParam(path + "Output/Param.csv");
centersX = csvread(path + "Output/CellCenters_X.csv");
centersP = csvread(path + "Output/CellCenters_P.csv");
time = stepNo * param.Dt;

% Generate titles
titleLine1 = "Plot of Solution " + solnName;
titleLine2 = "time = " + num2str(time) + "s";

% Graph all solution values, with ghost cell values
% fig = graphSoln(solnFilename, centersX, centersP, titleLine1, titleLine2);
% figure

% No ghost cells
centersX_noGhost = centersX(2:end-1, 2:end-1);
centersP_noGhost = centersP(2:end-1, 2:end-1);
if viewAngle == false
    fig = graphSoln(solnFilename, centersX_noGhost, centersP_noGhost, titleLine1, titleLine2);
else
    fig = graphSoln(solnFilename, centersX_noGhost, centersP_noGhost, titleLine1, titleLine2, viewAngle);
end

% Print pdf
if saveToPDF
    filename = "Output/" + solnName + "_" + num2str(time) + "s";
    saveas(fig, char(filename), 'pdf');
end

contourLevels = selectContourLevels(solnName);
figure
if viewAngle == false
    graphContour(solnFilename, centersX_noGhost, centersP_noGhost, contourLevels, titleLine1, titleLine2);
else
    graphContour(solnFilename, centersX_noGhost, centersP_noGhost, contourLevels, titleLine1, titleLine2, viewAngle);
end

