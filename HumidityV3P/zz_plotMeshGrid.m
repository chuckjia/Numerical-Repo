function zz_plotMeshGrid(meshGridX, meshGridP, cellCentersX, cellCentersP, plotCellCenters)
%PLOTMESHGRID Plot the mesh grid and cell centers in a 2D graph

if nargin < 5
    plotCellCenters = true;
end

mesh(meshGridX, meshGridP, zeros(size(meshGridX)));
hold on

% Graph the ghost cells
ghostCellColor = [218, 227, 243] ./ 255;  % Color used for the ghost cells. Light blue = [218, 227, 243]

% Shade all ghost cells along the four sides
[nrow, ncol] = size(meshGridX);
mesh(meshGridX(1:2,:), meshGridP(1:2,:), zeros(2, ncol), 'FaceColor', ghostCellColor);  % Left ghost cells
mesh(meshGridX((nrow-1):nrow,:), meshGridP((nrow-1):nrow,:), zeros(2, ncol), 'FaceColor', ghostCellColor);  % Right ghost cells
mesh(meshGridX(:,1:2), meshGridP(:,1:2), zeros(nrow, 2), 'FaceColor', ghostCellColor);  % Bottom ghost cells
mesh(meshGridX(:,(ncol-1):ncol), meshGridP(:,(ncol-1):ncol), zeros(nrow, 2), 'FaceColor', ghostCellColor);  % Top ghost cells

if plotCellCenters
    plot3(cellCentersX, cellCentersP, zeros(size(cellCentersP)), '.', 'MarkerSize', 10);
end
view(2)
hold off

end
