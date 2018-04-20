function mat = getVecFromFile( path )
% This function reads a 1D array from file

fileID = fopen(path);
mat = textscan(fileID, "%f"); 
fclose(fileID);

mat = cell2mat(mat);

end

