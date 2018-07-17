function levels = selectContourLevels(soln)
% SELECTCONTOURLEVELS

if soln == "T"
    levels = [255, 261, 273, 279, 284, 290, 293, 295, 298, 299, 305];
elseif soln == "q"
    levels = 1e-3 * [1, 3, 5, 8.83, 10.7, 12.5, 14, 16, 18, 20, 21];
else
    levels = 0;
end

end