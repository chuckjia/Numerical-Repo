function [cellCentersX, cellCentersP] = buildCellCenters(meshGridX, meshGridP)
%CALCCELLCENTERS Calculate the coordinates of the cell centers. 
%   For an illustration of the trapezoid cell, see he documentation file.
%
%   INPUT:: 
%      meshGridX: 2D matrix of size Nx+3 by Np+3. Its (i+1,j+1) entry is the x-coordinate x_{i-1/2}
%      meshGridP: 2D matrix of size Nx+3 by Np+3. Its (i+1,j+1) entry is the p-coordinate p_{i-1/2, j-1/2} 
% 
%   OUTPUT:: 
%      cellCentersX: 2D matrix of size Nx+2 by Np+2. Its (i+1,j+1) entry is the x-coordinate x_{i}
%      cellCentersP: 2D matrix of size Nx+2 by Np+2. Its (i+1,j+1) entry is the p-coordinate p_{i,j}
%

numCellsX = size(meshGridX, 1) - 1;  numCellsP = size(meshGridX, 2) - 1;
cellCentersX = zeros(numCellsX, numCellsP);
cellCentersP = zeros(numCellsX, numCellsP);

for i = 1:numCellsX
    for j = 1:numCellsP
        x_left = meshGridX(i, j);
        x_right = meshGridX(i + 1, j);
        p_bottLeft = meshGridP(i, j);
        p_bottRight = meshGridP(i + 1, j);
        p_topLeft = meshGridP(i, j + 1);
        p_topRight = meshGridP(i + 1, j + 1);
        
        x1 = (2 * x_left + x_right) / 3;
        p1 = (p_bottLeft + p_bottRight + p_topLeft) / 3;
        x2 = (x_left + 2 * x_right) / 3;
        p2 = (p_bottLeft + p_bottRight + p_topRight) / 3;        
        x3 = (x_left + 2 * x_right) / 3;
        p3 = (p_bottRight + p_topLeft + p_topRight) / 3;
        x4 = (2 * x_left + x_right) / 3;
        p4 = (p_bottLeft + p_topLeft + p_topRight) / 3;
        
        [cellCentersX(i, j), cellCentersP(i, j)] = calcLineIntersection(x1, p1, x2, p2, x3, p3, x4, p4);
    end
end

end

function [x, p] = calcLineIntersection(x1, p1, x2, p2, x3, p3, x4, p4) 
%CALCLINEINTERSECTION Calculate coordinate of the intersection point of two lines. 
%   Each line is represented by 2 points on that line. We assume that the points (x1, p1) and (x3, p3) are on 
%       line 1. The other 2 points are on line 2.
%

k1 = (p3 - p1) / (x3 - x1);
k2 = (p4 - p2) / (x4 - x2);
x = (k1 * x1 - p1 - k2 * x2 + p2) / (k1 - k2);
p = k1 * (x - x1) + p1;

end
