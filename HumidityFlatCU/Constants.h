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
#include "MathFcns.h"

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Numerical Scheme Constants
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

/*
 * User-set constants
 */
const int Nx = 100;  // Number of divisions on the x-direction
const int Np = 100;  // Number of divisions on the p-direction
const int numTimeSteps = 100;  // Number of time steps
const double Dt = 0.1;  // Size of time steps

/*
 * Auto-generated constants
 */
const int finalTime = numTimeSteps * Dt;

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Model Constants
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

/*
 * User-set constants
 */
const double x0 = 0;  // Leftmost x
const double xL = 50000;  // Rightmost x
const double pA = 200;
const double pB = 1000;

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

double u_fcn(double x, double p) {
	return 7.5 + cos(M_PI * p / p0_CONST) * cos(2 * M_PI * x / xL);
}

double omega_fcn(double x, double p) {
	return sin(M_PI * p / p0_CONST) * cos(2 * M_PI * x / xL);
}

double f_fcn(double input, double x, double p) {
	return u_fcn(x, p) * input;
}

double g_fcn(double input, double x, double p) {
	return omega_fcn(x, p) * input;
}

#endif /* CONSTANTS_H_ */
