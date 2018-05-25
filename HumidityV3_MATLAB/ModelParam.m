classdef ModelParam
    %MODELPARAM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        x0
        xf
        pA
        pB_x0
        pB_xf
        Nx
        Np
        Dt
        numTimeStep
        movieFrameFreq
        aveFreq
    end
    
    methods
        function obj = ModelParam(vec)
            if (length(vec) == 11)
                obj.x0 = vec(1);
                obj.xf = vec(2);
                obj.pA = vec(3);
                obj.pB_x0 = vec(4);
                obj.pB_xf = vec(5);
                obj.Nx = vec(6);
                obj.Np = vec(7);
                obj.Dt = vec(8);
                obj.numTimeStep = vec(9);
                obj.movieFrameFreq = vec(10);
                obj.aveFreq = vec(11);
            end
        end
    end
    
end

