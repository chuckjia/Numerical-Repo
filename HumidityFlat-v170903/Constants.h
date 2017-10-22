/*
 * Constants.h
 *
 *  Created on: Sep 3, 2017
 *      Author: chuckjia
 */

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#include <stdio.h>
#include <math.h>

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * General Settings
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Selection of finite volume method and model
const int fluxMethod = 0;  // 0: Godunov, 1: Classical FV
const int modelNum = 1;  // Model number

/* ----- ----- ----- ----- ----- -----
 * Scheme specifications
 * ----- ----- ----- ----- ----- ----- */

// Number of space divisions in both x and p directions in space
const int numDivisions = 300;
// Number of time steps
// Final time of numerical scheme
const double finalTime = 1e-3;
const int numTimeSteps = 100;
// Size of one time step
const double Dt = finalTime / numTimeSteps;

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Other Constants
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Geometry of domain
const double x0 = 0;  // x0 from model: x coordinate of left side of domain
const double xf = 100;  // xL from model: x coordinate of right side of domain
const double pA = 0;  // pA from model: p coordinate of bottom side of domain
const double pB = 100;  // pB from model: p coordinate of top side of domain

const int Nx = numDivisions;  // Number of divisions in x direction in space
const int Np = numDivisions;  // Number of divisions in p direction in space
const int numCellsX = Nx + 2;  // Number of cells in the x direction
const int numCellsP = Np + 2;  // Number of cells in the p direction

// Following constants are set up for code readability
const int lastRealIndexX = Nx;  // x direction index of cells on the right side of domain
const int lastRealIndexP = Np;  // p direction index of cells on the top side of domain
const int lastGhostIndexX = Nx + 1;  // x direction index of ghost cells on right side of domain
const int lastGhostIndexP = Np + 1;  // p direction index of ghost cells on top side of domain

const double Dx = (xf - x0) / Nx;  // x step size
const double Dp = (pB - pA) / Np;  // p step size
const double cellVol = Dx * Dp;  // Volume of a cell

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Getter Functions For the Mesh
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Return the x coordinate of the left side of cell (i, j)
double getCellLeftX(int i, int j) {
	return x0 + (i - 1) * Dx;
}

// Return the x coordinate of the right side of cell (i, j)
double getCellRightX(int i, int j) {
	return x0 + i * Dx;
}

// Return the p coordinate of the bottom side of cell (i, j)
double getCellBottP(int i, int j) {
	return pA + (j - 1) * Dp;
}

// Return the p coordinate of the top side of cell (i, j)
double getCellTopP(int i, int j) {
	return pA + j * Dp;
}

// Return the x coordinate of the center of cell (i, j)
double getCellCenterX(int i, int j) {
	return x0 + (i - 0.5) * Dx;
}

// Return the p coordinate of the center of cell (i, j)
double getCellCenterP(int i, int j) {
	return pA + (j - 0.5) * Dp;
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Model Functions
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Coefficients c1 and c2 used in the u function
const double c1_uFcn = M_PI / pB, c2_uFcn = 2 * M_PI / xf;

// Velocity function u
double u_fcn(double x, double p) {
	return 7.5 + cos(c1_uFcn * p) * cos(c2_uFcn * x);
}

// x derivative of the u function
double u_xDer_fcn(double x, double p) {
	return -cos(c1_uFcn * p) * c2_uFcn * sin(c2_uFcn * x);
}

// Coefficients c1 and c2 used in the w function
const double c1_wFcn = M_PI / pB, c2_wFcn = 2 * M_PI / xf;

// Velocity function w
double w_fcn(double x, double p) {
	return sin(c1_wFcn * p) * cos(c2_wFcn * x);
}

// p derivative of the w function
double w_pDer_fcn(double x, double p) {
	return c1_wFcn * cos(c1_wFcn * p) * cos(c2_wFcn * x);
}

#endif /* CONSTANTS_H_ */
