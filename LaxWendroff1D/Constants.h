/*
 * Constants.h
 *
 *  Created on: Oct 8, 2017
 *      Author: chuckjia
 */

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#include <stdio.h>
#include <math.h>

const int modelNum = 1;
const int methodNum = 0;  // 0 upwind Godonuv, 1 Lax-Wendroff, 2 Lax-Wendroff with flux limiter
const bool printMovie = false;
const bool printStat = false;

// Constants
const int numDivisions = 1000;
const double x0 = 0;
const double xf = 10;
const double a = 1;  // Velocity coefficient a

const double finalTime = 10;
const int numTimeSteps = 1000;
const double Dt = finalTime / numTimeSteps;

//const double Dt = 0.01;
//const int numTimeSteps = 800;
//const double finalTime = Dt * numTimeSteps;


// Constants from calculation
const double Dx = (xf - x0) / numDivisions;
const int numCells = numDivisions + 2;
const int lastRealIndex = numCells - 2;
const int lastGhostIndex = numCells - 1;

// Numerical solution
double u[numCells];

double getCellCenter(int i) {
	return x0 + (i - 0.5) * Dx;
}

#endif /* CONSTANTS_H_ */
