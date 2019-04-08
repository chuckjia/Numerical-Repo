function fig = graphContour(filename, centersX, centersP, contourLevels, titleLine1, titleLine2, viewAngle, proportion, colmapName)
%PLOTCONTOUR Summary of this function goes here
%   Detailed explanation goes here

if nargin < 9
    colmapName = 'parula';
    if nargin < 8
        proportion = 1;
    end
end

soln = csvread(filename);
graphSize = size(centersX);
solnSize = size(soln);

if isequal(graphSize + 2, solnSize)
    soln = soln(2:graphSize(1)+1, 2:graphSize(2)+1);
    graphSize = graphSize - 2;
elseif ~isequal(graphSize, solnSize)
    error("Mesh grid dimension size does not agree with solution matrix!");
    exit;
end

pRange = max(floor(graphSize(2) * (1-proportion)), 1) : graphSize(2);
[~, fig] = contourf(centersX(:, pRange), centersP(:, pRange), soln(:, pRange), contourLevels);
xlabel("x-axis");  ylabel("Solution Value");
caxis([contourLevels(1), contourLevels(end)])
colormap(colmapName)
colorbar

if nargin == 5
    title(titleLine1);
elseif nargin >= 6
    title({titleLine1, titleLine2, ""});
end

if nargin < 7
    viewAngle = false;
end

if length(viewAngle) == 2
    view(viewAngle);
end

end

