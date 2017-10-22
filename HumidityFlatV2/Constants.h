/*
 * Constants.h
 *
 *  Created on: Jul 10, 2017
 *      Author: chuckjia
 */

#ifndef CONSTANTS_H_
#define CONSTANTS_H_
#include "BasicFcns.h"

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * General Settings
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

/* ----- ----- ----- ----- ----- -----
 * Scheme specifications
 * ----- ----- ----- ----- ----- ----- */

// Number of space steps
const int numDivisions = 50;
// Number of time steps
const int numTimeSteps = 400;
// Size of time steps
const double finalTime = 1;  // Final time of numerical scheme
const double Dt = finalTime / numTimeSteps;

/* ----- ----- ----- ----- ----- -----
 * Scheme and model selection
 * ----- ----- ----- ----- ----- ----- */

// Scheme selection
const int numericalScheme = 0;  // 0 == Godunov, 1 == Central Upwind
const int timeScheme = 1;  // 1 == forward Euler, 2 == RK2
// Model selection
const int modelNumber = 1;  // 0 represents the original model

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Geometrical Constants
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

const double x0 = 0;  // Leftmost x
const double xf = 50000;  // Rightmost x: model value 5e4
const double pA = 200;  // Bottom p: model value 200
const double pB = 1000;  // Top p: model value 1000

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Auto-generated Variables
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

const int Nx = numDivisions;  // Number of divisions on the x-direction
const int Np = numDivisions;  // Number of divisions on the p-direction
const double Dx = (xf - x0) / Nx;  // Size of each x step
const double DxInv = Nx / (xf - x0);
const double Dp = (pB - pA) / Np;  // Size of each p step
const double DpInv = Np / (pB - pA);
const int numCellsX = Nx + 2;  // Number of cells in the x direction (including flat control volumes)
const int numCellsP = Np + 2;  // Number of cells in the p direction (including flat control volumes)
const int numSidesX = Nx + 1;  // The rightmost cell index on the x-direction (used for convenience)
const int numSidesP = Np + 1;  // The top cell index on the p-direction (used for convenience)
const int lastGhostX = Nx + 1;  // The top cell index on the p-direction (used for convenience)
const int lastGhostP = Np + 1;  // The top cell index on the p-direction (used for convenience)

#endif /* CONSTANTS_H_ */
