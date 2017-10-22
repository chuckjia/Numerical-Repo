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
const int fluxMethod = 1;  // 0: Godunov
const int modelNum = 1;  // Model number

/* ----- ----- ----- ----- ----- -----
 * Scheme specifications
 * ----- ----- ----- ----- ----- ----- */

// Number of space divisions in both x and p directions in space
const int numDivisions = 100;
// Number of time steps
const int numTimeSteps = 100;
// Size of one time step
const double Dt = 1e-5;
// Final time of numerical scheme
const double finalTime = numTimeSteps * Dt;

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Other Constants
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Geometry of domain
const double x0 = 0;  // x0 from model: x coordinate of left side of domain
const double xf = 10;  // xL from model: x coordinate of right side of domain
const double pA = 0;  // pA from model: p coordinate of bottom side of domain

// pB from model: the boundary on top side of domain
double pB(double x) {
	return 1000 - 200 * exp(-pow((x - 25000) / 3000, 2));
}

const int Nx = numDivisions;  // Number of divisions in x direction in space
const int Np = numDivisions;  // Number of divisions in p direction in space
const int numCellsX = Nx + 2;  // Number of cells in the x direction
const int numCellsP = Np + 2;  // Number of cells in the p direction

// Following constants are set up for code readability
const int lastRealIndexX = Nx;  // x direction index of cells on the right side of domain
const int lastRealIndexP = Np;  // p direction index of cells on the top side of domain
const int lastGhostIndexX = Nx + 1;  // x direction index of ghost cells on right side of domain
const int lastGhostIndexP = Np + 1;  // p direction index of ghost cells on top side of domain
const double numGridPtsX = numCellsX + 1;  // Number of grid points in the x direction (including ghosts)
const double numGridPtsP = numCellsP + 1;  // Number of grid points in the p direction (including ghosts)

const double Dx = (xf - x0) / Nx;  // x step size

#endif /* CONSTANTS_H_ */
