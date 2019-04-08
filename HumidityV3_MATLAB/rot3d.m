function rot3d(varargin)
%ROT3D Rotate the 3d plot as animation
%   This function would rotate the 3d plot around the z-axis. No input arguments would assume rotation uses default
%       parameters as specified below, with no videos saved.
%   
%   OPTIONAL INPUT ARGUMENTS:: 
%       'xyAngle':         A vector specifying the range of rotation angles. Default value: -135:1:225
%       'zAngle':          A scalar specifying the viewing angle from the z-axis. Default value: 30
%       'SaveVideo':       A logical value indicating whether video is recorded. Default value: false
%       'Filename':   A string specifying the filename of the saved video. Default value: "rotVideo.avi". The
%                              saved file will be stored in the folder named "Output/" within the project directory
%       'FrameRate':  A scalar specifying the video frame rate. Default value: 4
%

xyAngleVec = -135:1:225;
zAngle = 30;

% The default setting
if nargin == 0
    for az = xyAngleVec
        view(az, 30)
        drawnow
    end
    return
end

saveVideo = false;
Filename = "";
FrameRate = -1;
default_Filename = "rotVideo.avi";
default_FrameRate = 4;

p = inputParser;
validRotVec = @(x) isnumeric(x) && isvector(x);
addParameter(p, 'xyAngle', xyAngleVec, validRotVec);
validSingleNumber = @(x) isnumeric(x) && isscalar(x);
addParameter(p, 'zAngle', zAngle, validSingleNumber);
addParameter(p, 'SaveVideo', saveVideo, @islogical);
validFilename = @(x) isstring(x) || ischar(x);
addParameter(p, 'Filename', Filename, validFilename);
addParameter(p, 'FrameRate', FrameRate, validSingleNumber);
parse(p, varargin{:});

xyAngleVec = p.Results.xyAngle;
zAngle = p.Results.zAngle;
saveVideo = p.Results.SaveVideo;
Filename = p.Results.Filename;
if ischar(Filename)
    Filename = string(Filename);
end
FrameRate = p.Results.FrameRate;

if ~(Filename == "" && FrameRate < 0)  % If one of them is specified, then we will save video
    saveVideo = true;
end
if Filename == ""
    Filename = default_Filename;
end
if FrameRate < 0
    FrameRate = default_FrameRate;
end

% axis vis3d  % Can fix limit of z axis, but does not have good result

frameNo = 1;
for az = xyAngleVec
    view(az, zAngle)
    % zlim([0.005, 0.025])  % Use this to keep z axis fixed
    drawnow
    
    if saveVideo
        F(frameNo) = getframe(gcf);
        if frameNo == 1
            F = repmat(F, 1, length(xyAngleVec));
        end
        frameNo = frameNo + 1;
    end
end
% axis vis3d

if saveVideo
    movieName = "Output/" + Filename;
    v = VideoWriter(char(movieName));
    v.FrameRate = FrameRate;
    open(v);
    writeVideo(v, F);
    close(v);
    fprintf("Successfully saved video as " + movieName + "\n");
end

end

