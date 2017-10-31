function mat = getVecFromFile_fcn( folder, filename )

% This function reads a 1D array from file

repository = "/Users/chuckjia/Documents/Workspace/Git/Numerical-Repo/";
project = "HumidityV1/";

fileID = fopen(repository + project + folder + filename + ".txt");
mat = textscan(fileID, "%f"); 
fclose(fileID);

mat = cell2mat(mat);

end

