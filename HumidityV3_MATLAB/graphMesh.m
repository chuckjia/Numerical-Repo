function graphMesh( gridFile_X, gridFile_P, centersFile_X, centersFile_P )
%GRAPHMESH Summary of this function goes here
%   Detailed explanation goes here

plotCenter = false;
if (nargin == 4)
    plotCenter = true;
elseif (nargin ~= 2)
    error("Incorrect number of arguments");
end

% Read data from CSV files
gridX = csvread(gridFile_X);
gridP = csvread(gridFile_P);
assert(all( size(gridX) == size(gridP) ));
[nrow, ncol] = size(gridX);

% Graph the mesh
mesh(gridX, gridP, zeros(nrow, ncol));
hold on

% Graph the ghost cells
ghostCellColor = [218, 227, 243] ./ 255;  % Color used for the ghost cells. Light blue = [218, 227, 243]

% Shading all ghost cells along the four sides
mesh(gridX(1:2,:), gridP(1:2,:), zeros(2, ncol), 'FaceColor', ghostCellColor);  % Left ghost cells
mesh(gridX((nrow-1):nrow,:), gridP((nrow-1):nrow,:), zeros(2, ncol), 'FaceColor', ghostCellColor);  % Right ghost cells
mesh(gridX(:,1:2), gridP(:,1:2), zeros(nrow, 2), 'FaceColor', ghostCellColor);  % Bottom ghost cells
mesh(gridX(:,(ncol-1):ncol), gridP(:,(ncol-1):ncol), zeros(nrow, 2), 'FaceColor', ghostCellColor);  % Top ghost cells

% Graph the cell centers
if (plotCenter)
    centersX = csvread(centersFile_X);
    centersP = csvread(centersFile_P);
    assert(all( size(centersX) == size(centersP) ));
    assert(all( size(centersX) + 1 == size(gridX) ));
    
    plot3( centersX, centersP, zeros(size(centersX)), '.', 'MarkerSize', 10 );
end

view(2)
xlabel('x-axis');  ylabel('p-axis');
hold off

end

