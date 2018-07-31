# Humidity Project Documentation

This code performs numerical simulations for inviscid hydrostatic primitive equation models in the paper *Numerical Simulations of the Humid Atmosphere Above A Mountain*.

## Change Logs
#### Version 3.180524

- Changed calculation method of `Dp` in the function `calcCellCenterDp`, defined in `Mesh.h`. The calculation is now consistent with the domain mesh, instead of the actual physical domain. This affects the integration on p in the projection method. See code for details.


## How to Run the Code

To run the code, you can either compile the driver file or simply use the included makefile.

* The recommended way is to use makefile. In terminal, navigate to the project directory and use
  ~~~~
  make
  ~~~~
  This will compile and run the code.

* To manually compile, find the driver file and compile with C++ compilers. For example, if using `g++`, use the following to compile
  ~~~~
  g++ Driver.cpp -o output_filename
  ~~~~
  The compiled C++ file can be found under the directory `Output`.
  To execute the compiled binary file, use
  ~~~~
  ./output_filename 0
  ~~~~
  The argument `0` here suppresses printout of multi-line progress information. It does not affect computation results. The command `./output_filename` would work just as fine computation-wise, except that it might produce many lines of progress information. This is designed to work with Eclipse-style optimization on compilation of the code.


## Change Model Settings

In general, to change model settings, e.g. time step size, mesh size, etc., change variable values in the `Settings.h` file. To write new model equations, use the file `Models.h`.

In the file `Settings.h`:
* `numDivision` controls the size of the mesh on the physical domain. For example, if `numDivision` is set to be 50, then the model uses a mesh of size 50x50. Currently, we assume that the number of divisions on the x- and p- directions are the same. If needed, this can be easily changed by setting `Nx` and `Np` to different values in `Mesh.h`.

* `timeMethod` controls the time method. A value of 1 indicates that the first order forward Euler methods will be used. A value of 2 means the RK2 method is used. And a value of 4 means the RK4 method is used.

* `numTimeStep` is the number of time steps in the time advancement. `Dt` is the size of a single time step. And `finalTime` is the final time in the computation, calculated from the previous two variables.

* `aveSolnFreq` controls the frequency of averaging. For example, a value of 20 means the averaging method is applied to T every 20 time steps. The averaging method is applied to u, w, and phi every time step.

* `numProgMsg` controls the frequency of console printouts on current progress. For example, a value of 500 indicates during the computation, a total of 500 messages will be printed to console in the style of
~~~~
Current progress: 1.20 %, step no. 10
~~~~

* `movieFrameFreq` controls the frequency of writing solutions to files. Each of such files is a `.csv` file representing a numerical solution as a matrix. A value of 100 would indicate that the program writes to file every 100 time steps. These files will serve as the source of movie frames when making a movie of the simulation.
Note that a negative value indicates no such files are to be written during the life of the program.


* `calcL2NormFreq` controls the frequency of calculating and show L2 norms. A value of 100 indicates L2 norms is calculated and shown every 100 steps. A negative value means no L2 norms will be calculated during simulation.



In the file `models.h`:
* Two models are given. Model 0 is the physical case, while model 1 is the analytical case. All functions with suffix `_MDL0` are written for model 0, and functions with `_MDL1` are for model 1.

* Model and domain functions are constructed using function pointers to increase generality of code. To add new models, set function pointer values in the functions `setModels` and `setConditions`, which select model functions according to `modelNo` in `Settings.h`. The latter function, `setConditions`, is in the file `Conditions.h`. For example, the pB function, which describes the boundary geometry along the mountain surface, is represented by the function pointer `pB_fptr`. A function named `pB_fcn_MDL0` is written in `models.h` for the physical model as its pB function. Then when `modelNo` from `Settings.h` is set to be `0`, the value of `pB_fptr` will be set to `&pB_fcn_MDL0` in the function `setModels`.







## Computation Results

The results of the numerical simulation are stored in 2 directories: `Output` and `MovieFrames`.

* The following results are printed to `.csv` files in the directory `Output`:
    * the parameters used in the simulation: time and space step sizes, mesh size, etc.,
    * the L2-norms of numerical solutions during the time evolution, and
    * the numerical solution, the exact solution and the numerical errors at the final time.


* The numerical solutions during some of the time steps are printed to `.csv` files in the directory `MovieFrames`. These solutions can be used to generate surface and contour plots of the solutions, as well as movies, using the graphing code.



## Naming Conventions

#### Naming Conventions for Variables

##### In General
In general, variables have camel style names. For example, the following variable names are used in the project: `numTimeStep`, `qsVal`.

##### Constants
Constant variables have all-caps names or names with `_CONST` at the end. For example, the following constant variables are used in the project: `TWO_PI`, `g_CONST`.

##### Global Arrays
Global array variables have underscores at the end of their names. For example, `T_`, `meshGridP_`, etc.

##### Optional Variables
Variables serving as options have underscores before and after their names: `_[VAR_NAME]_`. For example, the variable `_aveResults_` is a `bool` variable that controls if averaging method is applied during simulations.

##### Miscellaneous
- Number of items: If a variable represents the number of certain items, then its name is usually shortened as `numTheItemName`. Here, we choose to use single forms of item names to save space. For example, `numTimeStep`, `numDivision`, etc. Most of the time, as a convention, we do not use `numTimeSteps` or `numDivisions`.




#### Naming Conventions for Functions

##### In General
In general, function names use camel styles. For example, `getCellTopRightP`, `runTimeSteps`.
They usually start with lowercase letters. However, this is not strictly followed. If it is more convenient to the start with a capital letter in the function name, then it is not recommendeded to convert it to uppercase. For example, the two functions `L_fcn` and `F_fcn` in the file `Constants.h` correspond to two physical functions named L and F. Then their natural names are used in the function naming.
