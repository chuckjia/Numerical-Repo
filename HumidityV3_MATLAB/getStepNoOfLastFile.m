function stepNo = getStepNoOfLastFile(filePath, solnName)
%GETSTEPNOOFLASTFILE Summary of this function goes here

solnName = char(solnName);
filename = getLastFileInDir(filePath, solnName);
stepNo = str2num(filename(3:end-4));

end

