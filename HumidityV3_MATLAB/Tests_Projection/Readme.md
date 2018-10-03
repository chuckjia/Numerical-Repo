# Tests On The Projection Method

This is a program designed to test the projection method. This program build a mesh for the spacial domain, assigns initial values to u on the mesh grid, and calculates the lambda_x values by solving a linear system as in (3.33) and (3.35) in the article [BCHTT15]. Then the projected u is calculated as u = initial_u - lambda_x.

## How to Run the Program

The program is coded in MATLAB.

- To run the program, simply run the `main(Nx)` function.

- The `main` function assumes a single input `Nx`, which represents number of cells in the x-direction. In this simulation, we assume that we use the same number of steps in the two spacial directions. That is, we have `Nx` = `Np`.

- The output of the function has two components `[u_afterProj, u_beforeProj]`, both matrices, representing the discrete values of u after and before projection.

- The console printout of the function gives the location x_i at which the derivative of\int_{pA}^{pB} u dp is the largest (in absolute value) and the value of the derivative at x_i.


## Some Extra Scripts and Functions

Two extra plotting scripts are provided.

- `plotMeshGrid` is a function that can plot the mesh grid. It provides a way to perform a sanity-check on the mesh grid. The function takes 4 matrices as input `meshGridX`, `meshGridP`, `cellCentersX`, and `cellCentersP`, representing the x- and p- coordinates of the mesh grids and the x- and p- coordinates of the cell centers, repectively.
