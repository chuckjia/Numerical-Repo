# Humidity Project Graphing Code Documentation

This MATLAB code makes plots using results computed from the HumidityV3 program.


## How to Run the Code

The main functions are `graphSoln`, `graphContour`, and `genMovie`.
Three examples are given on how to use these functions for graphics.

* The file `Eg_GraphSoln.m` gives an example of a graphing script.
  This script produces both the surface plot and the contour plot for the numerical solution at a certain time step.
  Before running the script, change the values of the variables `projectPath`, `stepNo`, and `solnName` accordingly.
  `projectPath` refers to the the path to the directory `HumidityV3`.
  `stepNo` is the step number for which we are generating plots.
  `solnName` is the name of the solution that we need to graph, i.e. `T`, `q`, `u`, or `w`.

* The file `Eg_GraphSoln_simple.m` gives two simpler examples on graphing numerical results. These examples only produce the surface plots.
Same as above, before running the script, change the values of the variables `projectPath`, `stepNo`, and `solnName` accordingly.

* The file `Eg_MakeMovie1.m` gives an example of generating a movie for the numerical simulation.
The variables `projectPath` and `solnName` have similar meanings as above.
The variable `steps` is a vector of the step numbers of all the steps which we want to include in the movie.

All the results printed to files can be found in the `Output` folder.
