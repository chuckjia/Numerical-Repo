/*
 * B_Constants.h
 *
 *  Created on: Apr 5, 2018
 *      Author: chuckjia
 */

#ifndef B_CONSTANTS_H_
#define B_CONSTANTS_H_

#include "A_Settings.h"

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Mathematical and Physical Constants
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Mathematical
double TWO_PI_CONST = 2 * M_PI;
double ONE_THIRD_CONST = 1. / 3.;
double ONE_SIXTH_CONST = 1. / 6.;

// Defined on page 100
double R_CONST = 287.;
double Rv_CONST = 461.5;
double Cp_CONST = 1004.;
double g_CONST = 9.8, gInv_CONST = 1. / 9.8;
double p0_CONST = 1000., p0Inv_CONST = 0.001;

// Defined on page 115: physical case
double T0_CONST = 300.;
double DeltaT_CONST = 50.;

// For convenience
double halfDt = 0.5 * Dt;
double oneSixthDt = Dt / 6.;

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Math Functions
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Returns the sign of a number
double sign_fcn(double x) { return x > 0 ? 1 : (x < 0 ? -1 : 0); }

// Empty function as place holders
void empty_fcn() { }

// Placeholder for a solution function
double zero_fcn(double x, double p, double t) {
	return 0;
}


#endif /* B_CONSTANTS_H_ */
