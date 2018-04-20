function F = genMovie( solnName, path, param, centersX, centersP, steps, saveToFile )
%GENMOVIE Summary of this function goes here
%   Detailed explanation goes here

contourLevels = selectContourLevels(solnName);
titleLine1 = "Plot of solution " + solnName;

frameNo = 0;
for stepNo = steps
    frameNo = frameNo + 1;
    time = stepNo * param.Dt;
    titleLine2 = "time = " + num2str(time) + "s";
    solnFilename = genSolnFilename(path, solnName, stepNo);
    
    viewAngle = [0, -90];
    graphContour(solnFilename, centersX, centersP, contourLevels, titleLine1, titleLine2, viewAngle);
    
    F(frameNo) = getframe(gcf);
    if (frameNo == 1) % Allocate memory space at the beginning
        F = repmat(F, 1, length(steps));
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

