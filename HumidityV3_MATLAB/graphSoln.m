function fig = graphSoln( filename, centersX, centersP, titleLine1, titleLine2, viewAngle )
%GRAPHSOLN Summary of this function goes here
%   Detailed explanation goes here

soln = csvread(filename);
graph_size = size(centersX);
soln_size = size(soln);

if isequal(graph_size + 2, soln_size)
    soln = soln(2:graph_size(1)+1, 2:graph_size(2)+1);
elseif ~isequal(graph_size, soln_size)
    error("Mesh grid dimension size does not agree with solution matrix!");
    exit;
end

fig = surf(centersX, centersP, soln);

if nargin == 4
    title(titleLine1);
elseif nargin >= 5
    title({titleLine1, titleLine2, ""});
    if nargin == 6
        view(viewAngle);
    end
end


end

