
%% 

clearAllScp
% projectPath = defaultProjectPath();  % Path to the outmost folder
projectPath = "~/Documents/Workspace/DataStorage/Humidity/HumidityV3/2018_11_29_15_39_25_good";
projectPath = genFolderPathName(projectPath);

% ===== ===== ===== =====
% Settings
% ===== ===== ===== ===== 

solnName = "q";
stepNo = -1;  % A value of -1 indicates graph the latest solution

resultFolder = "MovieFrames/";  % Commonly used: "MovieFrames/" or "Output/"
plotName = "Numerical Solution";  % "Solution", "Error", or "Exact Solution"

viewAngle = [0, -90];  % Viewing angle
contourProportion = 0.5;
colmapName = 'default'; % myGrayMap, 'default' or 'gray'

% Setting some of the values: Do Not Change
stepNo = genActualStepNo(stepNo, projectPath + resultFolder, solnName);

param = readParam(projectPath + "Output/Param.csv");
solnFileFullPath = genSolnFileFullPath(projectPath, resultFolder, solnName, stepNo);

% Mesh grid: cell centers
centersX = csvread(projectPath + "Output/CellCenters_X.csv");
centersP = csvread(projectPath + "Output/CellCenters_P.csv");
T_soln   = csvread(projectPath + "MovieFrames/T_" + num2str(stepNo) + ".csv");
q_soln   = csvread(projectPath + "MovieFrames/q_" + num2str(stepNo) + ".csv");
centersX = centersX(2:end-1, 2:end-1);  
centersP = centersP(2:end-1, 2:end-1); 
T_soln   = T_soln(2:end-1, 2:end-1);    
q_soln   = q_soln(2:end-1, 2:end-1);
qs_mat   = qs(T_soln, centersP);



% Generate graph titles
titleLine1 = "Plot of " + plotName + " " + solnName;
plotTime = param.Dt * stepNo;  displayPlotTime = plotTime;
titleLine2 = "Time = " + num2str(displayPlotTime) + "s";
contourLevels = selectContourLevels(solnName);

% q_soln = q_soln .* (q_soln <= qs_mat);

graphSize = size(centersX)-2;
solnSize = size(q_soln);
pRange = max(floor(graphSize(2) * (1-contourProportion)), 1) : graphSize(2);

contourf(centersX(:, pRange), centersP(:, pRange), q_soln(:, pRange), contourLevels);

xlabel("x-axis");  ylabel("Solution Value");
caxis([contourLevels(1), contourLevels(end)])
colormap(colmapName)
colorbar

title({titleLine1, titleLine2, ""});
view(viewAngle);


%% 

function val = qs(T, p)

val = 3.801664 .* exp(17.67 * (T - 273.15) ./ (T - 29.65)) ./ p;

end