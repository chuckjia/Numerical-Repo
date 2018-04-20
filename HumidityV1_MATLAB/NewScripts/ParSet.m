classdef ParSet
    %PARSET Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        x0
        xf
        pA
        pBx0
        pBxf
        pB
        Nx
        Np
        Dt
        numTimeSteps
        cMovieFrameFreq
        aveFreq
        
        Dx = (xf - x0) / Nx; Dp = (pB - pA) / Np;
        numCellsX = Nx + 2; numCellsP = Np + 2;
    end
    
    methods
    end
    
end

