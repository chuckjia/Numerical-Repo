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
using namespace std;

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Select Model and Scheme
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

int modelNo = 0;  // Model number: 0 = original physical model; other number = test models
int timeMethod = 4;  // Time method: 1 = Forward Euler, 2 = RK2, 4 = RK4


/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Scheme Specifications
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

const int numDivision = 200;  // Number of space divisions in both x and p directions

int numTimeStep = 100;  // Number of time steps
double Dt = 0.01;  // Size of one time step

double finalTime = numTimeStep * Dt;  // Final time


/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Settings
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Averaging method
int aveSolnFreq = 25;  // The frequency of using the averaging method. 0 or negative values indicate no averaging

// Movie I/O
int movieFrameFreq = 10000;  // The frequency of printing results to file as movie frames


bool _calcL2err_ = false;  // Choose whether to calculate, show, and print to file the L2 errors during computation

// Test cases
bool _printResultToFile_ = true;  // Choose if print numerical SOLUTION and ERRORS to file at the END of computation
bool _printExactSolnToFile_ = false;  // Choose if print EXACT solutions to file at the END of computation


#endif /* A_SETTINGS_H_ */
