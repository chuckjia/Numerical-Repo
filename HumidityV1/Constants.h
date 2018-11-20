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

int modelNo = 1;  // Model number. 0 = original physical model; other number = test models
int fluxMethod = 0;  // Select flux calculation method. 0 = Upwind Godunov

// Select time method.
// 1 = Forward Euler, 2 = RK2, 4 = RK4, 22 = BDF2
int timeMethod = 4;

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Scheme Specifications
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

const int numDivision = 100;  // Number of space divisions in both x and p directions

int numTimeStep = 100;  // Number of time steps
double Dt = 0.01;  // Size of one time step
double finalTime = numTimeStep * Dt;  // Final time

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Settings
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

/* ----- ----- ----- ----- ----- -----
 * Physical Simulations
 * ----- ----- ----- ----- ----- ----- */

// The frequency of using the average method
bool _aveResult_ = true;
int aveSolnFreq_T = 20;
// Choose whether to calculate, show, and print to file the L2 errors during computation
bool _calcL2err_ = true;

/* ----- ----- ----- ----- ----- -----
 * For testing cases
 * ----- ----- ----- ----- ----- -----*/

// Choose whether to print numerical SOLUTION and ERRORS to file at the END of computation
bool _printResultToFile_ = true;
// Choose if print EXACT solutions to file at the END of computation
bool _printExactSolnToFile_ = false;

/* ----- ----- ----- ----- ----- -----
 * Settings for movie I/O
 * ----- ----- ----- ----- ----- -----*/

// The frequency of printing results to file for movie frames
int movieFrameFreq = 1000;

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Set All Parameter in This File
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void setConstants() {
	// Empty for now
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Testing
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Print out messages for testing purposes
void tm() {
	printf("\n===== ===== ===== ===== ===== \n");
	printf("The program passed here.\n===== ===== ===== ===== =====\n");
}

#endif /* CONSTANTS_H_ */
