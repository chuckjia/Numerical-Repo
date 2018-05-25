function param = readParam( filename )
%READPARAM Summary of this function goes here
%   Detailed explanation goes here

vec = csvread(filename);
param = ModelParam(vec);

end

