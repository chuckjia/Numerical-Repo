
%% Test on the mesh grid: compare c and Matlab results
% Result 10/2: No discrepencies at all
clear; clc; close all

% Mesh and physical parameters
Nx = 200;  Np = Nx;  x0 = 0;  xf = 75000;  pA = 250; 

% From c
dataFolder = "test_data/";
gridX_c = csvread(dataFolder + "MeshGrid_X.csv");
gridP_c = csvread(dataFolder + "MeshGrid_P.csv");
gridX_c = gridX_c(2:end-1, 2:end-1);
gridP_c = gridP_c(2:end-1, 2:end-1);

% From Matlab
[gridX_m, gridP_m] = buildMesh(x0, xf, pA, Nx, Np);

% Discrepencies
discrepency_x = abs(gridX_c - gridX_m);
discrepency_p = abs(gridP_c - gridP_m);
max_discrepency_x = max(max(discrepency_x));
max_discrepency_p = max(max(discrepency_p));

small_num = 1e-15;
num_discrepency_x = sum(sum(discrepency_x > small_num));
num_discrepency_p = sum(sum(discrepency_p > small_num));

% Print results
fprintf("The max discrepencies in abs values:\n");
fprintf("  - in x is %1.6e,\n  - in p is %1.6e\n", max_discrepency_x, max_discrepency_p);
fprintf("The number of discrepencies:\n");
fprintf("  - in x is %d, %1.2f%%\n  - in p is %d, %1.2f%%.\n", ... 
    num_discrepency_x, num_discrepency_x / Nx^2 * 100, num_discrepency_p, num_discrepency_p / Nx^2 * 100);



%% Test on the mesh cell centers: compare c and Matlab results
clear; clc; close all

% Mesh and physical parameters
Nx = 200;  Np = Nx;  x0 = 0;  xf = 75000;  pA = 250; 

% From c
dataFolder = "test_data/";
centersX_c = csvread(dataFolder + "CellCenters_X.csv");
centersP_c = csvread(dataFolder + "CellCenters_P.csv");
centersX_c = centersX_c(2:end-1, 2:end-1);
centersP_c = centersP_c(2:end-1, 2:end-1);

% From Matlab
[gridX_m, gridP_m] = buildMesh(x0, xf, pA, Nx, Np);
[centersX_m, centersP_m] = calcCellCenters(gridX_m, gridP_m);

% Discrepencies
discrepency_x = abs(centersX_c - centersX_m);
discrepency_p = abs(centersP_c - centersP_m);
max_discrepency_x = max(max(discrepency_x));
max_discrepency_p = max(max(discrepency_p));

small_num = 1e-15;
num_discrepency_x = sum(sum(discrepency_x > small_num));
num_discrepency_p = sum(sum(discrepency_p > small_num));

% Print results
fprintf("The max discrepencies in abs values:\n");
fprintf("  - in x is %1.6e,\n  - in p is %1.6e\n", max_discrepency_x, max_discrepency_p);
fprintf("The number of discrepencies:\n");
fprintf("  - in x is %d, %1.2f%%\n  - in p is %d, %1.2f%%.\n", ... 
    num_discrepency_x, num_discrepency_x / Nx^2 * 100, num_discrepency_p, num_discrepency_p / Nx^2 * 100);
discrepency_x > small_num

%% Test on lambda values: compare c and Matlab results
clear; clc; close all

% Mesh and physical parameters
Nx = 200;  Np = Nx;  x0 = 0;  xf = 75000;  pA = 250; 

% From c
dataFolder = "test_data/";
lambdax_c = csvread(dataFolder + "lambdax.csv")';

% From matlab
[~, ~, ~, ~, ~, ~, ~, ~, lambdax_m] = main(Nx);
discrepencies = abs(lambdax_m - lambdax_c);
max(max(discrepencies))





%% Test on u values: compare c and Matlab results
clear; clc; close all

% Mesh and physical parameters
Nx = 200;  Np = Nx;  x0 = 0;  xf = 75000;  pA = 250; 

% From c
dataFolder = "test_data/";
u_before_c = csvread(dataFolder + "u_before.csv");
u_after_c = csvread(dataFolder + "u_after.csv");
u_before_c = u_before_c(2:end-1, 2:end-1);
u_after_c = u_after_c(2:end-1, 2:end-1);

% From Matlab
[u_after_m, u_before_m] = main(Nx); fprintf("\n");

% Discrepencies
discrepencies_before = abs(u_before_c - u_before_m);
discrepencies_after = abs(u_after_c - u_after_m);
u_before_max_discrepency = max(max(discrepencies_before));
u_after_max_discrepency = max(max(discrepencies_after));

small_num = 1e-15;
num_discrepencies_before = sum(sum(discrepencies_before > small_num));
num_discrepencies_after = sum(sum(discrepencies_after > small_num));

% Print results
fprintf("Before projection, the max discrepency in u is %1.8e.\n", u_before_max_discrepency);
fprintf("After projection, the max discrepency in u is %1.8e.\n", u_after_max_discrepency);
fprintf("Before projection, the number of discrepencies in u is %d.\n", num_discrepencies_before);
fprintf("After projection, the number of discrepencies in u is %d.\n", num_discrepencies_after);


%% Test on a,b,c coefficient values: compare c and Matlab results
clear; clc; close all

% Mesh and physical parameters
Nx = 200;  Np = Nx;  x0 = 0;  xf = 75000;  pA = 250; Dx = (xf - x0) / Nx;

% From c
dataFolder = "test_data/";
a_vec_c = csvread(dataFolder + "ab_vec_proj.csv");
b_vec_c = a_vec_c(2, :)';
a_vec_c = a_vec_c(1, :)';
c_vec_c = csvread(dataFolder + "c_vec_proj.csv")';
intU_vec_c = csvread(dataFolder + "intU_vec.csv")';

% From Matlab
[~,~,~,~, a_vec_m, b_vec_m, c_vec_m, intU_vec_m] = main(Nx);
a_vec_m = a_vec_m .* Dx - b_vec_m;
fprintf("\n");
discrepencies_a = abs(a_vec_c - a_vec_m);
[max_discrepency_a, max_loc_a] = max(discrepencies_a);
fprintf("Max discrepency of a is %1.5e, at location %d\n", max_discrepency_a, max_loc_a);
discrepencies_b = abs(b_vec_c - b_vec_m);
[max_discrepency_b, max_loc_b] = max(discrepencies_b);
fprintf("Max discrepency of b is %1.5e, at location %d\n", max_discrepency_b, max_loc_b);

discrepencies_c = abs(c_vec_c - c_vec_m .* Dx);
[max_discrepency_c, max_loc_c] = max(discrepencies_c);
fprintf("Max discrepency of c is %1.5e, at location %d\n", max_discrepency_c, max_loc_c);

discrepencies_intU = abs(intU_vec_c - intU_vec_m);
[max_discrepency_intU, max_loc_intU] = max(discrepencies_intU);
fprintf("Max discrepency of intU is %1.5e, at location %d\n", max_discrepency_intU, max_loc_intU);




%% 
clear; clc; close all

% Mesh and physical parameters
Nx = 200;  Np = Nx;  x0 = 0;  xf = 75000;  pA = 250; 

% From c
dataFolder = "test_data/";
centersX_c = csvread(dataFolder + "CellCenters_X.csv");
centersP_c = csvread(dataFolder + "CellCenters_P.csv");
centersX_c = centersX_c(2:end-1, 2:end-1);
centersP_c = centersP_c(2:end-1, 2:end-1);

% From Matlab
[gridX_m, gridP_m] = buildMesh(x0, xf, pA, Nx, Np);
[centersX_m, centersP_m] = calcCellCenters(gridX_m, gridP_m);

% Examine one cell
i = 13; 
j = 2;
x_grid_bottLeft = gridX_m(i, j);
x_grid_bottRight = gridX_m(i+1, j);
x_grid_topLeft = gridX_m(i, j+1);
x_grid_topRight = gridX_m(i+1, j+1);

p_grid_bottLeft = gridP_m(i, j);
p_grid_bottRight = gridP_m(i+1, j);
p_grid_topLeft = gridP_m(i, j+1);
p_grid_topRight = gridP_m(i+1, j+1);

fprintf("x_left = %1.5e, x_right = %1.5e\n", x_grid_bottLeft, x_grid_bottRight);
fprintf("p_bottLeft = %1.5e, p_bottRight = %1.5e, p_topRight = %1.5e, p_topLeft = %1.5e\n", ... 
    p_grid_bottLeft, p_grid_bottRight, p_grid_topRight, p_grid_topLeft);

x_center_m = centersX_m(i, j);
p_center_m = centersP_m(i, j);
fprintf("x_center = %1.30e, p_center = %1.30e\n", x_center_m, p_center_m);

% Write the result from c here
center_c = [4.68749999999999818101e+03,  2.55624999999999801048e+02];
discrepencies = [x_center_m, p_center_m] - center_c;
display(discrepencies)


%% Check compatibility conditions after projection in C++
clear; clc; close all

% Mesh and physical parameters
Nx = 200;  Np = Nx;  x0 = 0;  xf = 75000;  pA = 250;  Dx = (xf - x0) / Nx;

% From c
dataFolder = "test_data/";
u_after_c = csvread(dataFolder + "u_after.csv");
u_after_c = u_after_c(2:end-1, 2:end-1);

[meshGridX, meshGridP] = buildMesh(x0, xf, pA, Nx, Np);  % meshGridX/meshGridP are the x-/p-coordinates of the mesh grid
[cellCentersX, cellCentersP] = calcCellCenters(meshGridX, meshGridP);  % cellCentersX/cellCentersP are the x-/p-coordinates of the cell centers
Dp_vec = calcDp_cellCenter(x0, pA, Dx, Nx, Np, cellCentersX);

intU_vec_c = calcIntU(u_after_c, Dp_vec);
dxIntU_c = calcDxIntU(intU_vec_c, Dx);

x_vec = linspace(x0, xf, 200);
x_vec = x_vec(1:end-1);
% plot(x_vec, dxIntU_c)

hold on
[u_after_m, ~, ~, ~, ~, ~, ~, intU_vec_m] = main(Nx);
% plot(x_vec, dxIntU_c)
% hold off

intU_vec_m - intU_vec_c



%%
clear; clc
[u_afterProj, u_beforeProj, derAfterProj, derBeforeProj, ...
    a_vec, b_vec, c_vec, intU_vec, lambdax_vec] = main(200);

intU_vec



