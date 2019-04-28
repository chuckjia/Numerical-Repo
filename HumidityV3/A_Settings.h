/*
 * A_Settings.h
 *
 *  Created on: Oct 11, 2017
 *      Author: chuckjia
 *
 *  This file includes variables that controls the basic settings for the simulation.
 */

#ifndef A_SETTINGS_H_
#define A_SETTINGS_H_
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <string>
#include <cstdlib>
using namespace std;

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Select Model and Scheme
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

int modelNo    = 0;  // Model number: 0 = original physical model; other number = test models
int timeMethod = 4;  // Time method: 4 = RK4, 2 = RK2, 1 = Forward Euler, 0 = Control experiment with no time advancement


/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Scheme Specifications
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

/*
 * numDivision : Number of space divisions in both x and p directions
 * numTimeStep : Number of time steps
 * Dt          : Size of one time step
 *
 */

const int numDivision = 200;

int       numTimeStep = 40000;
double    Dt          = 0.5;

// Calculated
double    finalTime   = numTimeStep * Dt;  // Final time


/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Settings
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// aveSolnFreq    : The frequency of using the averaging method. 0 or negative values indicate no averaging
// numProgMsg     : Total number of progress messages. Progress messages are more frequent than movie frame prints and L2 norm info,
//                      guaranteed in setTimeSteps
// movieFrameFreq : Movie I/O: The frequency of printing results to file as movie frames
// calcL2NormFreq : Result evaluations: Choose whether to calculate, show, and print to file the L2 errors during computation
//
// movieFramesFolderName  : Name of the folder to print movie frames
// _printResultToFile_    : Choose if print numerical SOLUTION and ERRORS to file at the END of computation
// _printExactSolnToFile_ : Choose if print EXACT solutions to file at the END of computation
// _useCompatibleInitQ_   : Choose if use the more compatible initial conditions for q

int numMountain = 1;

int aveSolnFreq_T    = 32;  // Two mountains: 80, 3000, 2, 2, 2
int aveSolnFreq_q    = -1;
int aveSolnFreq_u    = 1;
int aveSolnFreq_w    = -1;
int aveSolnFreq_phix = -1;

int numProgMsg       = 200;
int movieFrameFreq   = 500;
int calcL2NormFreq   = 500;

string movieFramesFolderName = "MovieFrames";
bool _printResultToFile_    = true;
bool _printExactSolnToFile_ = false;
bool _useCompatibleInitQ_   = false;
bool _enforceTopBC_         = false;

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Validate Parameters
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void validateProgramParameters() {
	if (movieFrameFreq <= 0)    movieFrameFreq   = numTimeStep + 1;  // Do not print movie frames
	if (calcL2NormFreq <= 0)    calcL2NormFreq   = numTimeStep + 1;  // Do not calculate and show L2 errors

	if (aveSolnFreq_T <= 0)     aveSolnFreq_T    = numTimeStep + 1;
	if (aveSolnFreq_q <= 0)     aveSolnFreq_q    = numTimeStep + 1;
	if (aveSolnFreq_u <= 0)     aveSolnFreq_u    = numTimeStep + 1;
	if (aveSolnFreq_w <= 0)     aveSolnFreq_w    = numTimeStep + 1;
	if (aveSolnFreq_phix <= 0)  aveSolnFreq_phix = numTimeStep + 1;

	if (numMountain != 1 && numMountain != 2) {
		printf(">> ! Number of mountains is not set correctly. Using default value of 1 mountain.\n");
		numMountain = 1;
	}
	if (modelNo == 1) numMountain = 1;
}


#endif /* A_SETTINGS_H_ */
