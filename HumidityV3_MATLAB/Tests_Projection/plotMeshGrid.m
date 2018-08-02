function plotMeshGrid(meshGridX, meshGridP, cellCentersX, cellCentersP)
%PLOTMESHGRID Plot the mesh grid and cell centers in a 2D graph

mesh(meshGridX, meshGridP, zeros(size(meshGridX)));
hold on
plot3( cellCentersX, cellCentersP, zeros(size(cellCentersP)), '.', 'MarkerSize', 10 );
view(0, -90)
hold off

end

