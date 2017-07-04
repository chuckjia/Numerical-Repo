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

// Scheme selection
const int numericalScheme = 0;  // 0 represents Godunov, 1 represents Central Upwind
const int timeScheme = 2;  // 2 represents RK2, 4 represents RK4

// Space step size selection
const int numDivisions = 5;

// Time step size selection
const int numTimeSteps = 100;  // Number of time steps

// Test selection
const int testNumber = 10;  // Test 1, 2, or 3. 0 represents the original model

// Final Time
const double finalTime = 1;
const double Dt = finalTime / numTimeSteps;  // Size of time steps

/*const double Dt = 0.1;
const double finalTime = Dt * numTimeSteps;*/

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Numerical Scheme Constants
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

/*
 * Auto-generated constants
 */
const int Nx = numDivisions;  // Number of divisions on the x-direction
const int Np = numDivisions;  // Number of divisions on the p-direction

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Model Constants
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

/*
 * User-set constants
 */
const double x0 = 0;  // Leftmost x
const double xL = 5e4;  // Rightmost x: model value 5e4
const double pA = 2e2;  // Bottom p: model value 200
const double pB = 10e2;  // Top p: model value 1000

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

/*
 * Other (preset) model constants: constants used in the model functions
 */
const double theta_CONST = 1.25; // Parameter used in the van Leer's minmod limiters
const double p0_CONST = 1000;
const double R_CONST = 287;
const double Cp_CONST = 1004;
const double Rv_CONST = 461.50;

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Model Functions
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

/*
 * Coefficients for defining the function u (for the purpose of efficiency)
 */
const double u_fcn_COEFF1 = M_PI / p0_CONST,
		u_fcn_COEFF2 = 2 * M_PI / xL;
/*
 * The function u
 */
double u_fcn(double x, double p) {
	return -1;
	return 7.5 + cos(u_fcn_COEFF1 * p) * cos(u_fcn_COEFF2 * x);
}

/*
 * Coefficients for defining the function omega (for the purpose of efficiency)
 */
const double omega_fcn_COEFF1 = M_PI / p0_CONST,
		omega_fcn_COEFF2 = 2 * M_PI / xL;

/*
 * The function omega
 */
double omega_fcn(double x, double p) {
	return -1;
	return sin(omega_fcn_COEFF1 * p) * cos(omega_fcn_COEFF2 * x);
}

/*
 * The flux function f
 */
double f_fcn(double input, double x, double p) {
	return u_fcn(x, p) * input;
}

/*
 * The flux function g
 */
double g_fcn(double input, double x, double p) {
	return omega_fcn(x, p) * input;
}

#endif /* CONSTANTS_H_ */
