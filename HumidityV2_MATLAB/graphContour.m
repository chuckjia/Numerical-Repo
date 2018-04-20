function graphContour( filename, centersX, centersP, contourLevels, titleLine1, titleLine2, viewAngle )
%PLOTCONTOUR Summary of this function goes here
%   Detailed explanation goes here

soln = csvread(filename);
graph_size = size(centersX);
soln_size = size(soln);

if isequal(graph_size + 2, soln_size)
    soln = soln(2:graph_size(1)+1, 2:graph_size(2)+1);
    graph_size = graph_size - 2;
elseif ~isequal(graph_size, soln_size)
    error("Mesh grid dimension size does not agree with solution matrix!");
    exit;
end

p_range = floor(graph_size(2) * 0.7) : graph_size(2);
fig = contourf(centersX(:, p_range), centersP(:, p_range), soln(:, p_range), contourLevels);

if nargin == 5
    title(titleLine1);
elseif nargin >= 6
    title({titleLine1, titleLine2, ""});
    if nargin == 7
        view(viewAngle);
    end
end

end

