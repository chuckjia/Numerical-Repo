/*
 * Models.h
 *
 *  Created on: Oct 14, 2017
 *      Author: chuckjia
 */

#ifndef MODELS_H_
#define MODELS_H_
#include "Constants.h"

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Geometry of the Domain
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

double x0, xf, pA;
double (*pB_fcnPtr)(double x);  // Function pointer: pB function
double (*pBxDer_fcnPtr)(double x);  // Function pointer: derivative of pB function

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Mathematical and Physical Constants
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

double TWOPI_CONST = 2 * M_PI;
double R_CONST = 287;
double Rv_CONST = 461.50;
double Cp_CONST = 1004;
double g_CONST = 9.8, gInv_CONST = 1 / 9.8;

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Test Case 1: Test Case from the Section 4.1
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

/* ----- ----- ----- ----- ----- -----
 * Coefficients
 * ----- ----- ----- ----- ----- ----- */

// Coefficients used in the model
double xfCubed_coef_MDL1, xf6thPow_coef_MDL1;
double c1_pBxDer_coef_MDL1;
double pBx0_coef_MDL1;  // pB(x0)
double c1_exT_coef_MDL1;
double c1_exu_coef_MDL1;
double c1_exw_coef_MDL1;

/* ----- ----- ----- ----- ----- -----
 * Domain geometry
 * ----- ----- ----- ----- ----- ----- */

// pB function
double pB_fcn_MDL1(double x) {
	return 1000 - 200 * exp(-pow((x - 25000) / 3000, 2));
}

// The derivative of pB function
double pBxDer_fcn_MDL1(double x) {
	double a = (x - 25000) / 3000;
	return -c1_pBxDer_coef_MDL1 * a * exp(-a * a);
}

/* ----- ----- ----- ----- ----- -----
 * Manufacture solutions / IC
 * ----- ----- ----- ----- ----- ----- */

// Manufactured solution: exact T function
double exact_T_fcn_MDL1(double x, double p, double t) {
	return -p / R_CONST * cos(TWOPI_CONST * t) * x * pow(x - xf, 2) / xfCubed_coef_MDL1 * (
			c1_exT_coef_MDL1 * pow(p - (*pB_fcnPtr)(x), 2)
	);
}

// Manufactured solution: exact q function
double exact_q_fcn_MDL1(double x, double p, double t) {
	return 0;
}

// Manufactured solution: exact u function
double exact_u_fcn_MDL1(double x, double p, double t) {
	double a = p - pA, b = p - pBx0_coef_MDL1;
	return -c1_exu_coef_MDL1 * pow(x * (x - xf), 3) * (cos(TWOPI_CONST * t) + 20) *
			pow(a * b, 2) * (a + b);
}

// Manufactured solution: exact w function
double exact_w_fcn_MDL1(double x, double p, double t) {
	double a1 = x - xf,
			b1 = p - pB_fcn_MDL1(x); // b1 should be defined by pB_fcn_MDL1, NOT the pointer
	double a3 = pow(a1, 3), b3 = pow(b1, 3), x3 = pow(x, 3);
	return c1_exw_coef_MDL1 * (cos(TWOPI_CONST * t) + 20) * pow(p - pA, 3) * (
			x * x * a3 * b3
			+ a1 * a1 * b3 * x3
			- pBxDer_fcn_MDL1(x) * b1 * b1 * a3 * x3
	);
}

/* ----- ----- ----- ----- ----- -----
 * Source solutions
 * ----- ----- ----- ----- ----- ----- */

// Source function for the T equation
double source_T_fcn_MDL1(double T, double q, double u, double x, double p, double t) {
	return 0;
}

// Source function for the q equation
double source_q_fcn_MDL1(double T, double q, double u, double x, double p, double t) {
	return 0;
}

// Source function for the u equation
double source_u_fcn_MDL1(double T, double q, double u, double x, double p, double t) {
	return 0;
}

// Set all parameters in model 1
void setPar_MDL1() {
	// Parameters on the domain geometry
	x0 = 0;
	xf = 50000;
	pA = 200;
	pB_fcnPtr = &pB_fcn_MDL1;
	pBxDer_fcnPtr = &pBxDer_fcn_MDL1;

	// Function coefficients
	xfCubed_coef_MDL1 = pow(xf, 3);
	xf6thPow_coef_MDL1 = pow(xf, 6);
	c1_pBxDer_coef_MDL1 = 4. / 3.;
	pBx0_coef_MDL1 = (*pB_fcnPtr)(x0);
	c1_exT_coef_MDL1 = 2. / (450 * 450 * 450);
	c1_exu_coef_MDL1 = 3e-12 / xf6thPow_coef_MDL1;
	c1_exw_coef_MDL1 = 3e-12 / xf6thPow_coef_MDL1;
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Test Case 2: Constant Manufactured Solutions
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

/* ----- ----- ----- ----- ----- -----
 * Domain geometry
 * ----- ----- ----- ----- ----- ----- */

// pB function
double pB_fcn_MDL2(double x) {
	return 1000;
}

// Derivative of pB function
double pBxDer_fcn_MDL2(double x) {
	return 0;
}

/* ----- ----- ----- ----- ----- -----
 * Manufactured solutions / IC
 * ----- ----- ----- ----- ----- ----- */

// Exact T function
double exact_T_fcn_MDL2(double x, double p, double t) {
	return 1;
	// return sin(TWOPI_CONST * x / (xf - x0));
	// return sin(TWOPI_CONST * p / (pB_fcn_MDL2(x0) - pA));
}

// Exact q function
double exact_q_fcn_MDL2(double x, double p, double t) {
	return 1;
	// return sin(TWOPI_CONST * x / (xf - x0));
}

// Exact u function
double exact_u_fcn_MDL2(double x, double p, double t) {
	return 1;
	// return sin(TWOPI_CONST * x / (xf - x0));
}

// Exact w function
double exact_w_fcn_MDL2(double x, double p, double t) {
	return 0;
	// return sin(TWOPI_CONST * x / (xf - x0));
}

/* ----- ----- ----- ----- ----- -----
 * Source functions
 * ----- ----- ----- ----- ----- ----- */

double source_T_fcn_MDL2(double T, double q, double u, double x, double p, double t) {
	return 0;
}

double source_q_fcn_MDL2(double T, double q, double u, double x, double p, double t) {
	return 0;
}

double source_u_fcn_MDL2(double T, double q, double u, double x, double p, double t) {
	return 0;
}

void setPar_MDL2() {
	// Parameters on the geometry of the domain
	x0 = 0;
	xf = 50000;
	pA = 200;
	pB_fcnPtr = &pB_fcn_MDL2;
	pBxDer_fcnPtr = &pBxDer_fcn_MDL2;
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Set All Parameters
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Wrapper function to set all parameters for the selected model
void setModels() {
	// Default: if (modelNo == 1)
		setPar_MDL1();
	if (modelNo == 2)
		setPar_MDL2();
}

#endif /* MODELS_H_ */
