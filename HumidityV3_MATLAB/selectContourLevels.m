function levels = selectContourLevels(soln)
% SELECTCONTOURLEVELS

if soln == "T"
    n = 7;
    topLim = 299;
    bottLim = 280;
    levels = topLim - linspace(topLim - bottLim, 0, n); % Original
    % levels = [200, 220, 240, 260, 280, 300, 350, 400];
elseif soln == "q"
    n = 10;
    topLim = 0.0205;
    bottLim = 0.003;
    levels = topLim - linspace(topLim - bottLim, 0, n);
    % levels = 1e-3 * [1, 3, 5, 8.83, 10.7, 12.5, 14, 15.5, 18];
else
    levels = 0;
end

end