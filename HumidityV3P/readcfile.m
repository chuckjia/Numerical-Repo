function mat = readcfile(folderPath, filename)
%READCFILE Summary of this function goes here
%   Detailed explanation goes here

mat = csvread(fullfile(folderPath, folderPath));

end

