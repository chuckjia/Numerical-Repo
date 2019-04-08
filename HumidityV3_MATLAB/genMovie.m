function F = genMovie( solnName, path, param, centersX, centersP, steps, saveToFile )
%GENMOVIE Generate movies for the numerical results from the HumidityV3 program

contourLevels = selectContourLevels(solnName);
titleLine1 = "Plot of solution " + solnName;

figPositionVec = [0, 600, 900, 500];  % [left bottom width height]
set(gcf, 'PositionMode', 'Manual')

frameNo = 0;
skip = true;  % Bug on figure size temporary fix
for stepNo = [steps(1), steps]
    frameNo = frameNo + 1;
    time = stepNo * param.Dt;
    titleLine2 = "time = " + num2str(time) + "s";
    solnFilename = genSolnFilename(path, solnName, stepNo);
    
    viewAngle = [0, -90];
    proportion = 0.5;
    graphContour(solnFilename, centersX, centersP, contourLevels, titleLine1, titleLine2, viewAngle, proportion);
    set(gcf, 'Position', figPositionVec)
    F(frameNo) = getframe(gcf);
    
    if (frameNo == 1) % Allocate memory space at the beginning
        F = repmat(F, 1, length(steps));
    end
    
    if (frameNo == 2 && skip)  % Bug on figure size temporary fix
        frameNo = 1;
        skip = false;
        continue
    end
end

if nargin == 7 && saveToFile
    movieName = "Output/" + solnName + ".avi";
    v = VideoWriter(char(movieName));
    v.FrameRate = 4;
    open(v);
    writeVideo(v, F);
    close(v);
end

fprintf("Movie generation completed.\n");

end

