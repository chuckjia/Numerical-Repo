function [fig, soln] = graphSoln(filename, centersX, centersP, titleLine1, titleLine2, viewAngle)
%GRAPHSOLN Graph numerical solutions
%   Detailed explanation goes here

soln = csvread(filename);
graphSize = size(centersX);
solnSize = size(soln);

if isequal(graphSize + 2, solnSize)
    soln = soln(2:graphSize(1)+1, 2:graphSize(2)+1);
elseif ~isequal(graphSize, solnSize)
    error("Mesh grid dimension size does not agree with solution matrix!");
    exit;
end

fig = surf(centersX, centersP, soln);
xlabel('x-axis'); ylabel('p-axis');

if nargin == 4
    title(titleLine1);
elseif nargin >= 5
    title({titleLine1, titleLine2, ""});
end

if nargin < 6
    viewAngle = false;
end

if length(viewAngle) == 2
    view(viewAngle);
end

end

