function levels = selectContourLevels( soln )
%SELECTCONTOURLEVELS Summary of this function goes here
%   Detailed explanation goes here

if soln == "T"
    levels = [190, 267, 273, 279, 284, 290, 293, 295, 297, 297.8];
elseif soln == "q"
    levels = 1e-3 * [1, 3, 5, 8.83, 10.7, 12.5, 14, 15, 16];
else
    levels = 0;

end

