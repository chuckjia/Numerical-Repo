function F2 = slowMotion( F, rate, movieName )
%SLOWMOTION Summary of this function goes here
%   Detailed explanation goes here

numFrame = length(F) * rate;

F2 = repmat(F(1), 1, numFrame);

for i = 1:numFrame
    F2(i) = F(ceil(i / rate));
end


v = VideoWriter(movieName);
v.FrameRate = 15;
open(v);
writeVideo(v, F2);
close(v);

end

