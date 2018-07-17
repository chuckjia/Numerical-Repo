# Humidity Project Documentation

This code performs numerical simulations for inviscid hydrostatic primitive equation models in the paper *Numerical Simulations of the Humid Atmosphere Above A Mountain*.

## Change Logs
#### Version 3.180524

- Changed calculation method of `Dp` in the function `calcCellCenterDp`, defined in `Mesh.h`. The calculation is now consistent with the domain mesh, instead of the actual physical domain. This affects the integration on p in the projection method. See code for details.


## How to Run the Code

To run the code, you can either compile the driver file or simply use the included makefile.
* To manually compile, find the driver file and compile with C++ compilers. For example, if using `g++`, use the following to compile
  ~~~~
  g++ Driver.cpp -o output_filename
  ~~~~
  The compiled C++ file can be found under the directory `Output`.
  To execute the compiled binary file, use
  ~~~~
  ./output_filename
  ~~~~
  The name of the driver cpp file might vary slightly across versions.

* To use makefile, in terminal, navigate to the project directory and use
  ~~~~
  make
  ~~~~
  This will compile and run the code.

To change model settings, e.g. time step size, mesh size, etc., change variable values in the `Settings.h` file. To write new model equations, use the file `Models.h`.



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
