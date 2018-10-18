function F2 = slowMotion( F, rate, movieName )
%SLOWMOTION Summary of this function goes here
%   Detailed explanation goes here

numFrame = length(F) * rate;

F2 = repmat(F(1), 1, numFrame);

fprintf("Processing slow motion movie...\n");
for i = 1:numFrame
    F2(i) = F(ceil(i / rate));
end

fprintf("Saving slow motion movie...\n");
v = VideoWriter(movieName);
v.FrameRate = 15;
open(v);
writeVideo(v, F2);
close(v);

fprintf("Slow motion movie generated and saved.\n");

end

