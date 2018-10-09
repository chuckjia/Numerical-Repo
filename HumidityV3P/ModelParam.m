classdef ModelParam
    %MODELPARAM A class that holds all computation parameters
    %
    %   Author: Chuck Jia
    %   Date: 10/08/2018
    %
    
    properties
        % Domain geometry
        x0
        xf
        pA
        Dx
        
        % Mesh size
        Nx
        Np
        Dt
        Nt
        t_end
        
        % Other parameters
        AveRate
        MovieFrameRate
        NumMsg
    end
    
    methods
        % Class constructor
        function obj = ModelParam(x0, xf, pA, Nx, Np, Dt, Nt, AveRate, MovieFrameRate, NumMsg)
            % Validate parameters
            if Nx <= 0 || Np <= 0 || Dt <= 0 || Nt <= 0
                error("Nx, Np, Dt, and Nt should all be positive numbers!\n");
            end
            if MovieFrameRate > 0 && (Nt / NumMsg) > MovieFrameRate
                NumMsg = ceil(Nt / MovieFrameRate);
            end
            
            obj.x0 = x0;
            obj.xf = xf;
            obj.pA = pA;
            
            obj.Nx = Nx;
            obj.Np = Np;
            obj.Dt = Dt;
            obj.Nt = Nt;
            
            obj.MovieFrameRate = MovieFrameRate;
            obj.AveRate = AveRate;
            obj.NumMsg = NumMsg;
            
            obj.Dx = (xf - x0) / Nx;
            obj.t_end = Dt * Nt;
        end
        
        % Extract parameters
        function [x0, xf, pA, Nx, Np, Dx, Dt, Nt, t_end, AveRate, MovieFrameRate, NumMsg] = extractParam(obj)
            x0 = obj.x0;
            xf = obj.xf;
            pA = obj.pA;
            
            Nx = obj.Nx;
            Np = obj.Np;
            Dx = obj.Dx;
            Dt = obj.Dt;
            Nt = obj.Nt;
            t_end = obj.t_end;
            
            MovieFrameRate = obj.MovieFrameRate;
            AveRate = obj.AveRate;
            NumMsg = obj.NumMsg;
        end
        
        % Display contents
        function showParam(obj)
            fprintf("Model Parameters:\n");
            fprintf("  - Domain Geometry:   x0 = %1.2f,  xf = %1.2f,  pA = %1.2f\n", obj.x0, obj.xf, obj.pA);
            fprintf("  - Mesh Paramters:    Nx = %1.0f,  Np = %1.0f,  Dt = %1.2f,  Nt = %1.0f,  t_end = %1.2fs\n", ...
                obj.Nx, obj.Np, obj.Dt, obj.Nt, obj.t_end);
            
            fprintf("  - Averaging Rate:    ")
            if obj.AveRate > 0
                fprintf("every %1.0f steps\n", obj.AveRate);
            else
                fprintf("no averaging is applied\n")
            end
            
            fprintf("  - Movie Frame Rate:  ")
            if obj.MovieFrameRate > 0
                fprintf("every %1.0f steps\n", obj.MovieFrameRate)
            else
                fprintf("not recording\n")
            end
            
            fprintf("  - Progress Message:  ")
            if obj.MovieFrameRate > 0
                fprintf("every %1.0f steps\n", obj.NumMsg)
            else
                fprintf("all progress messages suppressed\n")
            end
            
            fprintf("\n");
        end
        
    end
    
end
