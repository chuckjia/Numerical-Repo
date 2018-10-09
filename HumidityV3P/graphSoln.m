function fig = graphSoln(soln, cellCentersX, cellCentersP, titleLine1, titleLine2, viewAngle)
%GRAPHSOLN Graph numerical solutions

graphSize = size(cellCentersX);
solnSize = size(soln);

if ~isequal(graphSize, solnSize)
    error("Mesh grid dimension size does not agree with solution matrix!");
    exit;
end

fig = surf(cellCentersX, cellCentersP, soln);
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

