/*
 * Constants.h
 *
 *  Created on: Jun 20, 2017
 *      Author: chuckjia
 *
 *  This file contains constants used in the model and the numerical schemes.
 */

#ifndef CONSTANTS_H_
#define CONSTANTS_H_
#include "BasicFcns.h"

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * General Settings
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Step size selection
const int numDivisions = 200;  // Number of space steps
const int numTimeSteps = 100;  // Number of time steps

// Test selection
const int testNumber = 1;  // Test 1, 2, or 3. 0 represents the original model

// Final Time

const double finalTime = 1;
const double Dt = finalTime / numTimeSteps;  // Size of time steps

/*const double Dt = 0.01;
const double finalTime = Dt * numTimeSteps;*/

// Scheme selection
const int numericalScheme = 0;  // 0 represents Godunov, 1 represents Central Upwind
const int timeScheme = 2;  // 1 represents forward Euler, 2 represents RK2

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Numerical Scheme Constants
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

/*
 * Auto-generated constants
 */
const int Nx = numDivisions;  // Number of divisions on the x-direction
const int Np = numDivisions;  // Number of divisions on the p-direction

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Mesh Constants
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

/*
 * User-set constants
 */
const double x0 = 0;  // Leftmost x
const double xL = 50000;  // Rightmost x: model value 5e4
const double pA = 200;  // Bottom p: model value 200
const double pB = 1000;  // Top p: model value 1000

/*
 * Auto-generated constants
 */
const double Dx = (xL - x0) / Nx;  // Size of each x step
const double DxInv = Nx / (xL - x0);
const double Dp = (pB - pA) / Np;  // Size of each p step
const double DpInv = Np / (pB - pA);
const int numCellsX = Nx + 2;  // Number of cells in the x direction (including flat control volumes)
const int numCellsP = Np + 2;  // Number of cells in the p direction (including flat control volumes)
const int lastIndexX = Nx + 1;  // The rightmost cell index on the x-direction (used for convenience)
const int lastIndexP = Np + 1;  // The top cell index on the p-direction (used for convenience)

#endif /* CONSTANTS_H_ */
