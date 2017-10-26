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
int modelNo = 2;  // Model number
int fluxMethod = 0;  // 0: Godunov

/* ----- ----- ----- ----- ----- -----
 * Scheme specifications
 * ----- ----- ----- ----- ----- ----- */

// Number of space divisions in both x and p directions in space
const int numDivisions = 100;
// Number of time steps
int numTimeSteps = 100;
// Size of one time step
double Dt = 1e-5;
// Final time of numerical scheme
double finalTime = numTimeSteps * Dt;

void msg() {
	printf("===== ===== ===== ===== ===== \n");
	printf("The program passed here.\n===== ===== ===== ===== =====\n");
}

#endif /* CONSTANTS_H_ */
