function param = showParam(projectPath)
%SHOWPARAM Summary of this function goes here
%   Detailed explanation goes here

if nargin == 0
    compAccID = whichcomputer();
    if compAccID == 1
        projectPath = "~/Documents/Workspace/git/Numerical-Repo/HumidityV3/";
    elseif compAccID == 2
        projectPath = "~/git/Numerical-Repo/HumidityV3/";
    else
        error("Unknown computer account. Please provide the project path.\n");
    end
end
    
param = readParam(projectPath + "Output/Param.csv");
    
fprintf("Mesh size = %dx%d, numTimeSteps = %d, Dt = %1.4f\n", ... 
    param.Nx, param.Np, param.numTimeStep, param.Dt);
end

