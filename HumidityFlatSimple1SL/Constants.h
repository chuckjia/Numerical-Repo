/*
 * Constants.h
 *
 *  Created on: Jul 27, 2017
 *      Author: chuckjia
 */

#ifndef CONSTANTS_H_
#define CONSTANTS_H_
#include <stdio.h>
#include <math.h>

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * General Settings
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

/* ----- ----- ----- ----- ----- -----
 * Scheme specifications
 * ----- ----- ----- ----- ----- ----- */

// Number of space steps
const int numDivisions = 100;
// Number of time steps
const int numTimeSteps = 100;
// Size of time steps
const double Dt = 1e-5;
const double finalTime = numTimeSteps * Dt;

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Other Constants
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

const double x0 = 0;
const double xf = 10;
const double pA = 0;
const double pB = 10;

const int Nx = numDivisions;
const int Np = numDivisions;
const int numCellsX = Nx + 2, numCellsP = Np + 2;
const int lastRealIndexX = Nx, lastRealIndexP = Np;
const int lastGhostIndexX = Nx + 1, lastGhostIndexP = Np + 1;

const double Dx = (xf - x0) / Nx, Dp = (pB - pA) / Np;
const double cellVol = Dx * Dp;

#endif /* CONSTANTS_H_ */
