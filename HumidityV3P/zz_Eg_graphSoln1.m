%% Examples to graph solutions
% This is file of examples on how to use the graphSoln function to plot the solutions

%% Example 1
% Simple example on graphing solutions

clc; close all

% ===== ===== ===== =====
% Settings
% ===== ===== ===== =====

solnName = "T";
timeStep = -1;
plotName = "Numerical Solution";  % "Solution", "Error", or "Exact Solution"

saveToPDF = false;  % Whether to save result to PDF files
viewAngle = [0, -90];  % Viewing angle
graphGhostCells = false;  % Whether to graph ghost cells

% Generate mesh
param = setModelParam();
[meshGridX, meshGridP] = buildMesh(param);
[cellCentersX, cellCentersP] = buildCellCenters(meshGridX, meshGridP);

% Validate parameters
if ischar(solnName)
    solnName = string(solnName);
end

if timeStep < 0
    timeStep = param.Nt;
end

switch solnName
    case "T"
        soln = T_;
    case "q"
        soln = q_;
    case "u"
        soln = u_;
    case "w"
        soln = w_;
    otherwise
        error("Solution does not exist.\n");
end

% Generate graph
titleLine1 = plotName + " " + solnName;
titleLine2 = "t = " + num2str(param.Dt * timeStep) + "s";

if graphGhostCells
    fig = graphSoln(soln, cellCentersX, cellCentersP, titleLine1, titleLine2, viewAngle);
else
    fig = graphSoln(soln(2:end-1,2:end-1), cellCentersX(2:end-1,2:end-1), cellCentersP(2:end-1,2:end-1), titleLine1, titleLine2, viewAngle);
end




