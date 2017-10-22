/*
 * Constants.h
 *
 *  Created on: Oct 11, 2017
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
int fluxMethod = 1;  // 0: Godunov
int modelNo = 1;  // Model number

/* ----- ----- ----- ----- ----- -----
 * Scheme specifications
 * ----- ----- ----- ----- ----- ----- */

// Number of space divisions in both x and p directions in space
const int numDivisions = 10;
// Number of time steps
int numTimeSteps = 10;
// Size of one time step
double Dt = 1e-5;
// Final time of numerical scheme
double finalTime = numTimeSteps * Dt;

#endif /* CONSTANTS_H_ */
