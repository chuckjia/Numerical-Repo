/*
 * Constants.h
 *
 *  Created on: Jun 10, 2017
 *      Author: chuckjia
 */

#ifndef CONSTANTS_H_
#define CONSTANTS_H_
#include "MathFcns.h"
/*
 * The following inclusions are used in other files. However,
 * they are included here for clarity
 */
#include "TestFcns.h"

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Model constants
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

const int Nx = 100;
const int Np = 100;
const int numTimeSteps = 100;
const double Dt = 0.1;
const int finalTime = numTimeSteps * Dt;

const double x0 = 0;
const double xL = 50;
const double Dx = (xL - x0) / Nx;
const int numGridPtsXDir = Nx + 1;
const int numGridPtsPDir = Np + 1;
const int numCellsXDir = Nx + 2;
const int numCellsPDir = Np + 2;

// Initialize the solution
double soln[numCellsXDir][numCellsPDir][2];

const double pA = 200;

double pB(double x) {
	return 1000. - 200. * exp(- pow(x - 25000., 2) / 9000000.);
}

double Dp(double x) {
	return (pB(x) - pA) / Np;
}

// Theta is the parameter used in the van Leer's minmod limiters
const double theta_CONST = 1.25;

// Constants used in the model functions

const double p0_CONST = 1000;
const double R_CONST = 287;
const double Cp_CONST = 1004;
const double Rv_CONST = 461.50;


/* ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Model functions
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

double uFcn(double x, double p) {
	return 7.5 + cos(M_PI * p / p0_CONST) * cos(2 * M_PI * x / xL);
}

double omegaFcn(double x, double p) {
	return sin(M_PI * p / p0_CONST) * cos(2 * M_PI * x / xL);
}

double fFcn(double indVar, double x, double p) {
	return uFcn(x, p) * indVar;
}

double gFcn(double indVar, double x, double p) {
	return omegaFcn(x, p) * indVar;
}

#endif /* CONSTANTS_H_ */
