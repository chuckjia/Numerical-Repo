function filePath = genSolnFileFullPath(path, folder, solnName, stepNo)
%GENSOLNFILEFULLPATH Summary of this function goes here
%   Detailed explanation goes here

filePath = path + folder + solnName + "_" + num2str(stepNo) + ".csv";

end

