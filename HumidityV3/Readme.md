# Humidity Project Documentation

This code performs numerical simulations for inviscid hydrostatic primitive equation models in the paper *Numerical Simulations of the Humid Atmosphere Above A Mountain*.

## Change Logs
#### Version 3.180524

- Changed calculation method of `Dp` in the function `calcCellCenterDp`, defined in `Mesh.h`. The calculation is now consistent with the mesh cells. See code for details.


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
  The driver cpp file might vary slightly across versions.

* To use makefile, in terminal, navigate to the project directory and use
  ~~~~
  make
  ~~~~
  This will compile and run the code.

To change model settings, e.g. time step size, mesh size, etc., change variable values in the `Settings.h` file. To write new model equations, use the file `Models.h`.





## Naming Conventions

#### Naming Conventions for Variables

##### In General
In general, variables have camel style names. For example, the following variable names are used in the project: `numTimeStep`, `qsVal`.

##### Constants
Constant variables have all-caps names or names with `_CONST` at the end. For example, the following constant variables are used in the project: `TWO_PI`, `g_CONST`.

##### Optional Variables
Variables serving as options have underscores before and after their names: `_[VAR_NAME]_`. For example, the variable `_aveResults_` is a `bool` variable that controls if averaging method is applied during simulations.

##### Miscellaneous
- Number of items: If a variable represents the number of certain items, then its name is usually shortened as `numTheItemName`. Here, we choose to use single forms of item names to save space. For example, `numTimeStep`, `numDivision`, etc. Most of the time, as a convention, we do not use `numTimeSteps` or `numDivisions`.




#### Naming Conventions for Functions

#### In General
In general, function names use camel styles. For example, `getCellTopRightP`, `runTimeSteps`.
