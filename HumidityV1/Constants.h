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
#include <time.h>
#include <iostream>
using namespace std;

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Select Model and Scheme
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

int modelNo = 1;  // Model number
int fluxMethod = 0;  // 0: Godunov

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Scheme Specifications
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

const int numDivisions = 250;  // Number of space divisions in both x and p directions
int numTimeSteps = 1000;  // Number of time steps
double Dt = 1e-7;  // Size of one time step
double finalTime = numTimeSteps * Dt;  // Final time

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Testing
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Print out messages for testing purposes
void testmsg() {
	printf("===== ===== ===== ===== ===== \n");
	printf("The program passed here.\n===== ===== ===== ===== =====\n");
}

#endif /* CONSTANTS_H_ */
