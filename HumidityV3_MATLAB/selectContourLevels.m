function levels = selectContourLevels(soln)
% SELECTCONTOURLEVELS

if soln == "T"
    levels = [255, 261, 273, 279, 285, 291, 294, 297, 300]; % Original
    % levels = [200, 220, 240, 260, 280, 300, 350, 400];
elseif soln == "q"
    levels = 1e-3 * [1, 3, 5, 8.83, 10.7, 12.5, 14, 15.5, 18, 20];
else
    levels = 0;
end

end