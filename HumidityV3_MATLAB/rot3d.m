function rot3d(varargin)
%ROT3D Rotate the 3d plot as animation
%   This function would rotate the 3d plot around the z-axis. No input arguments would assume rotation uses default
%       parameters as specified below, with no videos saved.
%   
%   OPTIONAL INPUT ARGUMENTS:: 
%       'xyAngle':         A vector specifying the range of rotation angles. Default value: -135:1:225
%       'zAngle':          A scalar specifying the viewing angle from the z-axis. Default value: 30
%       'SaveVideo':       A logical value indicating whether video is recorded. Default value: false
%       'VideoFilename':   A string specifying the filename of the saved video. Default value: "rotVideo.avi". The
%                              saved file will be stored in the folder named "Output/" within the project directory
%       'VideoFrameRate':  A scalar specifying the video frame rate. Default value: 4
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
videoFilename = "rotVideo.avi";
videoFrameRate = 4;

p = inputParser;
validRotVec = @(x) isnumeric(x) && isvector(x);
addParameter(p, 'xyAngle', xyAngleVec, validRotVec);
validSingleNumber = @(x) isnumeric(x) && isscalar(x);
addParameter(p, 'zAngle', zAngle, validSingleNumber);
addParameter(p, 'SaveVideo', saveVideo, @islogical);
addParameter(p, 'VideoFilename', videoFilename, @isstring);
addParameter(p, 'VideoFrameRate', videoFrameRate, validSingleNumber);
parse(p, varargin{:});

xyAngleVec = p.Results.xyAngle;
zAngle = p.Results.zAngle;
saveVideo = p.Results.SaveVideo;
videoFilename = p.Results.VideoFilename;
videoFrameRate = p.Results.VideoFrameRate;

frameNo = 1;
for az = xyAngleVec
    view(az, zAngle)
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
    movieName = "Output/" + videoFilename;
    v = VideoWriter(char(movieName));
    v.FrameRate = videoFrameRate;
    open(v);
    writeVideo(v, F);
    close(v);
end

end

