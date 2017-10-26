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

// Geometry of domain
double x0;  // x0 from model: x coordinate of left side of domain
double xf;  // xL from model: x coordinate of right side of domain
double pA;  // pA from model: p coordinate of bottom side of domain
double (*pB_fcnPtr)(double x);  // Function pointer for the pB function
double (*pBxDer_fcnPtr)(double x);  // Function pointer for the pB_x function

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Common Constants
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Common constants used in the model, which are physical or mathematical constants
double TWOPI_CONST = 2 * M_PI;
double R_CONST = 287;
double Rv_CONST = 461.50;
double Cp_CONST = 1004;
double g_CONST = 9.8, gInv_CONST = 1 / 9.8;

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Test Case 1
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

double xfCubed_coef_MDL1, xf6thPow_coef_MDL1;
double c1_pBxDer_coef_MDL1;
double pBx0_coef_MDL1;  // pB(x0)
double c1_exT_coef_MDL1;
double c1_exu_coef_MDL1;
double c1_exw_coef_MDL1;

// pB: the boundary on top side of domain
double pB_fcn_MDL1(double x) {
	return 1000 - 200 * exp(-pow((x - 25000) / 3000, 2));
}

// The derivative of pB
double pBxDer_fcn_MDL1(double x) {
	double a = (x - 25000) / 3000;
	return -c1_pBxDer_coef_MDL1 * a * exp(-a * a);
}

// Exact T function
double exact_T_fcn_MDL1(double x, double p, double t) {
	return -p / R_CONST * cos(TWOPI_CONST * t) * x * pow(x - xf, 2) / xfCubed_coef_MDL1 * (
			c1_exT_coef_MDL1 * pow(p - (*pB_fcnPtr)(x), 2)
	);
}

// Exact q function
double exact_q_fcn_MDL1(double x, double p, double t) {
	return 0;
}

// Exact u function
double exact_u_fcn_MDL1(double x, double p, double t) {
	double a = p - pA, b = p - pBx0_coef_MDL1;
	return -c1_exu_coef_MDL1 * pow(x * (x - xf), 3) * (cos(TWOPI_CONST * t) + 20) *
			pow(a * b, 2) * (a + b);
}

// Exact w function
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

double source_T_fcn_MDL1(double T, double q, double u, double x, double p, double t) {
	return 0;
}

double source_q_fcn_MDL1(double T, double q, double u, double x, double p, double t) {
	return 0;
}

double source_u_fcn_MDL1(double T, double q, double u, double x, double p, double t) {
	return 0;
}

void setPar_MDL1() {
	// Parameters on the geometry of the domain
	x0 = 0;
	xf = 50000;
	pA = 200;
	pB_fcnPtr = &pB_fcn_MDL1;
	pBxDer_fcnPtr = &pBxDer_fcn_MDL1;

	// Other coefficients in computations
	xfCubed_coef_MDL1 = pow(xf, 3);
	xf6thPow_coef_MDL1 = pow(xf, 6);
	c1_pBxDer_coef_MDL1 = 4. / 3.;
	pBx0_coef_MDL1 = (*pB_fcnPtr)(x0);
	c1_exT_coef_MDL1 = 2. / (450 * 450 * 450);
	c1_exu_coef_MDL1 = 3e-12 / xf6thPow_coef_MDL1;
	c1_exw_coef_MDL1 = 3e-12 / xf6thPow_coef_MDL1;
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Test Case 2
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// pB: the boundary on top side of domain
double pB_fcn_MDL2(double x) {
	return 1000;
}

// The derivative of pB
double pBxDer_fcn_MDL2(double x) {
	return 0;
}

// Exact T function
double exact_T_fcn_MDL2(double x, double p, double t) {
	return 1;
}

// Exact q function
double exact_q_fcn_MDL2(double x, double p, double t) {
	return 0;
}

// Exact u function
double exact_u_fcn_MDL2(double x, double p, double t) {
	return 1;
}

// Exact w function
double exact_w_fcn_MDL2(double x, double p, double t) {
	return 1;
}

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
 * Wrapper Function to Set Parameters
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// This is the wrapper function used to set parameters for the model.
void setModels() {
	// Default: if (modelNo == 1)
	setPar_MDL1();
	if (modelNo == 2)
		setPar_MDL2();
}

#endif /* MODELS_H_ */
