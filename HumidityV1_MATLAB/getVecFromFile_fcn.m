function mat = getVecFromFile_fcn( folder, filename )

% This function reads a 1D array from file

fileID = fopen(folder + filename + ".txt");
mat = textscan(fileID, "%f"); 
fclose(fileID);

mat = cell2mat(mat);

end

