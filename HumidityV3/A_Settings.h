/*
 * Constants.h
 *
 *  Created on: Oct 11, 2017
 *      Author: chuckjia
 */

#ifndef A_SETTINGS_H_
#define A_SETTINGS_H_
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <string>
using namespace std;

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Select Model and Scheme
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

int modelNo = 0;  // Model number. 0 = original physical model; other number = test models

// Select time method
// 1 = Forward Euler, 2 = RK2, 4 = RK4
int timeMethod = 4;

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Scheme Specifications
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

const int numDivisions = 200;  // Number of space divisions in both x and p directions

int numTimeStep = 20000;  // Number of time steps
double Dt = 0.5;  // Size of one time step
double finalTime = numTimeStep * Dt;  // Final time

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Settings
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

/* ----- ----- ----- ----- ----- -----
 * Physical Simulations
 * ----- ----- ----- ----- ----- ----- */

// The frequency of using the average method
bool averageResult_opt = true;
int aveFreq = 25;
// Choose whether to calculate, show, and print to file the L2 errors during computation
bool calcL2err_opt = false;

/* ----- ----- ----- ----- ----- -----
 * For testing cases
 * ----- ----- ----- ----- ----- ----- */

// Choose whether to print numerical SOLUTION and ERRORS to file at the END of computation
bool printResToFile_opt = true;
// Choose if print EXACT solutions to file at the END of computation
bool printExactSolnToFile_opt = false;

/* ----- ----- ----- ----- ----- -----
 * Settings for movie I/O
 * ----- ----- ----- ----- ----- -----*/

// The frequency of printing results to file for movie frames
int movieFrameFreq = 100;

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Set All Parameter in This File
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void setAllSettings() {
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

#endif /* A_SETTINGS_H_ */
