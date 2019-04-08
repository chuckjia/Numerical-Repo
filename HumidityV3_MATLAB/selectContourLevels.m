function levels = selectContourLevels(soln)
% SELECTCONTOURLEVELS

if soln == "T"
    
    
    
    % Original problem: [274, 299] for color and gray
    %     levels = [200, 220, 240, 260, 280, 300, 350, 400];
    
    % Two-mountain case
    %     n = 7;
    %     topLim = 299;
    %     bottLim = 274;
    %     levels = topLim - linspace(topLim - bottLim, 0, n);
    
    % Other test cases
    n = 7;
    topLim = 299;
    bottLim = 270;
    levels = topLim - linspace(topLim - bottLim, 0, n);
    
elseif soln == "q"
    
    % Original problem : [0.003, 0.0205] for color and gray
%         n = 8;
%         topLim = 0.02;
%         bottLim = 0.003;
%         levels = topLim - linspace(topLim - bottLim, 0, n);
    
    % Two-mountain case, n = 9
    %     n = 9;
    %     topLim = 0.0205;
    %     bottLim = -0.001;
    %     levels = topLim - linspace(topLim - bottLim, 0, n);
    
    % Other test cases
    n = 9;
    topLim = 0.01;
    bottLim = 0.0015;
    levels = topLim - linspace(topLim - bottLim, 0, n);
    
else
    levels = 0;
end

end