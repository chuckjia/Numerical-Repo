function filename = genSolnFilename( path, solnName, stepNo )
%GENSOLNFILENAME Summary of this function goes here
%   Detailed explanation goes here

if nargin == 3
    filename = path + solnName + "_" + int2str(stepNo) + ".csv";
elseif nargin == 2
    filename = path + solnName + ".csv";
end

end

