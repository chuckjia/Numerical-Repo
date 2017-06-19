/*
 * Constants.h
 *
 *  Created on: Jun 17, 2017
 *      Author: chuckjia
 */

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#include <stdio.h>
#include <math.h>
/*
 * The following inclusions are used in other files. However,
 * they are included here for clarity
 */
#include "MathFcns.h"
#include "TestFcns.h"

const int Nx = 10;
const int Np = 10;

const double x0 = 0;
const double xL = 50;
const double Dx = (xL - x0) / Nx;
const int numXGridPts = Nx + 1;
const int numPGridPts = Np + 1;
const int numCellsXDir = Nx + 2;
const int numCellsPDir = Np + 2;

// Initialize the solution
double soln[numCellsXDir][numCellsPDir][2];

const double pA = 200 * 0.1;

double pB(double x) {
	return (1000. - 200. * exp(- pow(x - 25000., 2) / 9000000.)) * 0.1;
}

double Dp(double x) {
	return (pB(x) - pA) / Np;
}

// Theta is the parameter used in the van Leer's minmod limiters
const double theta = 1.25;

const double p0 = 1000;

double uFcn(double x, double p) {
	return 7.5 + cos(M_PI * p / p0) * cos(2 * M_PI * x / xL);
}

double omegaFcn(double x, double p) {
	return sin(M_PI * p / p0) * cos(2 * M_PI * x / xL);
}

double fFcn(double indVar, double x, double p) {
	return uFcn(x, p) * indVar;
}

double gFcn(double indVar, double x, double p) {
	return omegaFcn(x, p) * indVar;
}

#endif /* CONSTANTS_H_ */
