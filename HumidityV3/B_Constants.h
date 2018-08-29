/*
 * B_Constants.h
 *
 *  Created on: Apr 5, 2018
 *      Author: chuckjia
 *
 *  This file contains constants and common functions used in mathematical and physical settings. Also included are common functions used
 *  in testing.
 */

#ifndef B_CONSTANTS_H_
#define B_CONSTANTS_H_

#include "A_Settings.h"

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Mathematical and Physical Constants
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Mathematical constants
double TWO_PI = 2 * M_PI;
double ONE_THIRD = 1. / 3.;
double ONE_SIXTH = 1. / 6.;

// Physical constants, defined on page 100
double R_CONST = 287.;
double Rv_CONST = 461.50;
double Cp_CONST = 1004.;
double g_CONST = 9.8;
double p0_CONST = 1000., p0Inv_CONST = 0.001;

// Physical constants used in the physical case, defined on page 115
double T0_CONST = 300.;
double DeltaT_CONST = 50.;

// For convenience
double halfDt = 0.5 * Dt;
double oneSixthDt = Dt / 6.;


/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Predefined Functions
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

/**
 * Mathematical Functions
 */

// Returns the sign of a number
double sgn_fcn(double x) { return x > 0 ? 1 : (x < 0 ? -1 : 0); }
int sgn_fcn(int x) { return x > 0 ? 1 : (x < 0 ? -1 : 0); }

// Placeholder for a solution function
double zero_fcn(double x, double p, double t) { return 0; }

/**
 * Physical Functions
 */

// Saturation specific humidity function. Used in the initial condition for q
double qs_fcn(double T, double p) {
	return 3.801664 / p * exp(17.67 * (T - 273.15) / (T - 29.65));  // Leading coefficient is 0.622 * 6.112
}

// Delta function defined on p.100, line 10-15
double delta_fcn(double q, double w, double qsVal) {
	return 0.25 * (1 - sgn_fcn(w)) * (1 + sgn_fcn(q - qsVal));
}

// L function defined on p.100, line 15-20
double L_fcn(double T) { // @suppress("Name convention for function")
	return 2.5008e6 - 2.3e3 * (T - 275.);
}

// F function defined in (2.2) on p.100
double F_fcn(double T, double qsVal, double LVal) { // @suppress("Name convention for function")
	return qsVal * T * (LVal * R_CONST - Cp_CONST * Rv_CONST * T) / (Cp_CONST * Rv_CONST * T * T + qsVal * LVal * LVal);
}

/**
 * Empty Placeholder Functions
 */

void empty_fcn() { }
void empty_fcn(int x) { }  // Used by aveSoln_fptr
void empty_fcn(double x) { }

/**
 * Testing Functions
 */

void printTimeUsed(clock_t start, clock_t end, string msg) {
	double timeUsed = ((double) (end - start)) / CLOCKS_PER_SEC;
	printf("- %s Time used = %1.2fs.\n", msg.c_str(), timeUsed);
}

void printTimeUsed(clock_t start, clock_t end, string unit, string msg) {
	double timeUsed = ((double) (end - start)) / CLOCKS_PER_SEC;
	if (unit == "s")
		printf("- %s Time used = %1.2fs.\n", msg.c_str(), timeUsed);
	else if (unit == "ms")
		printf("- %s Time used = %1.2fms.\n", msg.c_str(), timeUsed * 1000);
	else
		printf("Error: Wrong unit || In printTimeUsed with message %s\n", msg.c_str());
}

// Print out messages for testing purposes
void tttt() { printf("\n===== ===== ===== ===== ===== \n>> The program passed here.\n"); printf("===== ===== ===== ===== =====\n"); }


#endif /* B_CONSTANTS_H_ */
