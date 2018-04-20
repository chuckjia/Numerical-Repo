function cellCenters = getCellCentersFromFile( file_path, param )
%GETCELLCENTERSFROMFILE Summary of this function goes here
%   Detailed explanation goes here

assert(isa(param, 'ParSet'))

cellCenters = reshape(getVecFromFile(file_path), [param.Np, param.Nx]);

end

