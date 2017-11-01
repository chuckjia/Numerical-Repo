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
double R_CONST = 287.;
double Rv_CONST = 461.5;
double Cp_CONST = 1004.;
double g_CONST = 9.8, gInv_CONST = 1 / 9.8;
double T0_CONST = 300.;
double p0_CONST = 1000.;
double DeltaT_CONST = 50.;

double halfDt = 0.5 * Dt;

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
double c1_exT_coef_MDL1, c2_exT_coef_MDL1, c3_exT_coef_MDL1, c4_exT_coef_MDL1;
double c1_exu_coef_MDL1, c1_exu_pDer_coef_MDL1;
double c1_exw_coef_MDL1;

/* ----- ----- ----- ----- ----- -----
 * Domain geometry
 * ----- ----- ----- ----- ----- ----- */

// pB function
double pB_fcn_MDL1(double x) {
	return 1000 - 200 * exp(-pow((x - 25000) / 3000, 2));
}

// The derivative of pB function
double pB_xDer_fcn_MDL1(double x) {
	double tmp = (x - 25000) / 3000;
	return c1_pBxDer_coef_MDL1 * tmp * exp(-tmp * tmp);
}

/* ----- ----- ----- ----- ----- -----
 * Manufacture solutions / IC
 * ----- ----- ----- ----- ----- ----- */

// This is the exact T function excluding the t terms
double exact_T_helper_fcn_MDL1(double x, double p) {
	return c1_exT_coef_MDL1 * x * pow(x - xf, 2) * (
			p * (c2_exT_coef_MDL1 * pow(p - pB_fcn_MDL1(x), 2) - c3_exT_coef_MDL1)
			+ c4_exT_coef_MDL1
	);
}

// Manufactured solution: exact T function
double exact_T_fcn_MDL1(double x, double p, double t) {
	return cos(TWOPI_CONST * t) * exact_T_helper_fcn_MDL1(x, p);
}

// x-derivative of the exact T function
double exact_T_xDer_fcn_MDL1(double x, double p, double t) {
	double x_minus_xf = x - xf, pB_minus_p = pB_fcn_MDL1(x) - p;
	double xPart1 = x * x_minus_xf * x_minus_xf,
			xPart2 = p * (c2_exT_coef_MDL1 * pB_minus_p * pB_minus_p - c3_exT_coef_MDL1)
			+ c4_exT_coef_MDL1;
	return c1_exT_coef_MDL1 * cos(TWOPI_CONST * t) * (
			(3 * x - xf) * x_minus_xf * xPart2 +
			xPart1 * (2 * c2_exT_coef_MDL1 * pB_minus_p * pB_xDer_fcn_MDL1(x) * p)
	);
}

// The p-derivative of the exact T function
double exact_T_pDer_fcn_MDL1(double x, double p, double t) {
	double p_minus_pB = p - pB_fcn_MDL1(x);
	return c1_exT_coef_MDL1 * cos(TWOPI_CONST * t) * x * pow(x - xf, 2) * (
			c2_exT_coef_MDL1 * p_minus_pB * p_minus_pB - c3_exT_coef_MDL1 +
			p * c2_exT_coef_MDL1 * 2 * p_minus_pB
	);
}

// The t-derivative of the exact T function
double exact_T_tDer_fcn_MDL1(double x, double p, double t) {
	return -TWOPI_CONST * sin(TWOPI_CONST * t) * exact_T_helper_fcn_MDL1(x, p);
}

// Manufactured solution: exact q function
double exact_q_fcn_MDL1(double x, double p, double t) {
	return 0;
}

// The terms in u that involves only x and p
double exact_u_fcn_helper_MDL1(double x, double p) {
	double p_minus_pA = p - pA, p_minus_pB = p - pB_fcn_MDL1(x);
	return c1_exu_coef_MDL1 * pow(x * (x - xf), 3) *
			pow(p_minus_pA * p_minus_pB, 2) * (p_minus_pA + p_minus_pB);
}

// Manufactured solution: exact u function
double exact_u_fcn_MDL1(double x, double p, double t) {
	return (cos(TWOPI_CONST * t) + 20) * exact_u_fcn_helper_MDL1(x, p);
}

// The x-derivative of the exact u function
double exact_u_xDer_fcn_MDL1(double x, double p, double t) {
	double x_minus_xf = x - xf,
			p_minus_pA = p - pA, pB_minus_p = pB_fcn_MDL1(x) - p,
			pB_xDer_val = pB_xDer_fcn_MDL1(x);
	double xPart1 = pow(x, 3), xPart2 = pow(x_minus_xf, 3),
			xPart3 = pB_minus_p * pB_minus_p, xPart4 = p_minus_pA - pB_minus_p;
	return c1_exu_coef_MDL1 * p_minus_pA * p_minus_pA * (cos(TWOPI_CONST * t) + 20) * (
			3 * x * x * xPart2 * xPart3 * xPart4
			+ xPart1 * 3 * x_minus_xf * x_minus_xf * xPart3 * xPart4
			+ xPart1 * xPart2 * 2 * pB_minus_p * pB_xDer_val * xPart4
			- xPart1 * xPart2 * xPart3 * pB_xDer_val
	);
}

// The p-derivative of the exact u function
double exact_u_pDer_fcn_MDL1(double x, double p, double t) {
	double x_minus_xf = x - xf,
			p_minus_pA = p - pA, p_minus_pB = p - pB_fcn_MDL1(x);
	double p_sum_part = p_minus_pA + p_minus_pB, p_prod_part = p_minus_pA * p_minus_pB;
	return c1_exu_pDer_coef_MDL1 * pow(x * x_minus_xf, 3) * (cos(TWOPI_CONST * t) + 20)
			* p_prod_part * (p_sum_part * p_sum_part + p_prod_part);
}

// The t-derivative of the exact u function
double exact_u_tDer_fcn_MDL1(double x, double p, double t) {
	return -TWOPI_CONST * sin(TWOPI_CONST * t) * exact_u_fcn_helper_MDL1(x, p);
}

// Manufactured solution: exact w function
double exact_w_fcn_MDL1(double x, double p, double t) {
	double xPart1 = x * (x - xf), p_minus_pB = p - pB_fcn_MDL1(x);
	return c1_exw_coef_MDL1 * pow(p - pA, 3) * (cos(TWOPI_CONST * t) + 20)
			* xPart1 * xPart1 * p_minus_pB * p_minus_pB * (
					p_minus_pB * (2 * x - xf) - pB_xDer_fcn_MDL1(x) * xPart1
			);
}

// The p-derivative of the exact w function
double exact_w_pDer_fcn_MDL1(double x, double p, double t) {
	double xTerm = x * (x - xf), p_minus_pB = p - pB_fcn_MDL1(x), p_minus_pA = p - pA,
			two_x_minus_xf = 2 * x - xf;
	double pPart1 = pow(p_minus_pA, 3), pPart2 = p_minus_pB * p_minus_pB,
			pPart3 = p_minus_pB * two_x_minus_xf - pB_xDer_fcn_MDL1(x) * xTerm;
	return c1_exw_coef_MDL1 * xTerm * xTerm * (cos(TWOPI_CONST * t) + 20) * (
			3 * p_minus_pA * p_minus_pA * pPart2 * pPart3
			+ pPart1 * 2 * p_minus_pB * pPart3
			+ pPart1 * pPart2 * two_x_minus_xf
	);
}

/* ----- ----- ----- ----- ----- -----
 * Source solutions
 * ----- ----- ----- ----- ----- ----- */

// Source function for the T equation
double source_T_fcn_MDL1(double T, double q, double u, double x, double p, double t) {
	return exact_T_tDer_fcn_MDL1(x, p, t)
			+ exact_u_fcn_MDL1(x, p, t) * exact_T_xDer_fcn_MDL1(x, p, t)
			+ exact_w_fcn_MDL1(x, p, t) * exact_T_pDer_fcn_MDL1(x, p, t);
}

// Source function for the q equation
double source_q_fcn_MDL1(double T, double q, double u, double x, double p, double t) {
	return 0;
}

// Source function for the u equation
double source_u_fcn_MDL1(double T, double q, double u, double x, double p, double t) {
	return exact_u_tDer_fcn_MDL1(x, p, t)
			+ exact_u_fcn_MDL1(x, p, t) * exact_u_xDer_fcn_MDL1(x, p, t)
			+ exact_w_fcn_MDL1(x, p, t) * exact_u_pDer_fcn_MDL1(x, p, t);
}

// Set all parameters in model 1
void setPar_MDL1() {
	// Parameters on the domain geometry
	x0 = 0;
	xf = 50000;
	pA = 200;
	pB_fcnPtr = &pB_fcn_MDL1;
	pBxDer_fcnPtr = &pB_xDer_fcn_MDL1;

	// Function coefficients
	xfCubed_coef_MDL1 = pow(xf, 3);
	xf6thPow_coef_MDL1 = pow(xf, 6);
	c1_pBxDer_coef_MDL1 = 4. / 30.;
	pBx0_coef_MDL1 = (*pB_fcnPtr)(x0);

	c1_exT_coef_MDL1 = -1 / (g_CONST * xfCubed_coef_MDL1);
	c2_exT_coef_MDL1 = 3 * g_CONST / (450 * 450 * 450 * R_CONST);
	c3_exT_coef_MDL1 = DeltaT_CONST / p0_CONST;
	c4_exT_coef_MDL1 = DeltaT_CONST - T0_CONST;

	c1_exu_coef_MDL1 = -3e-12 / xf6thPow_coef_MDL1;
	c1_exu_pDer_coef_MDL1 = 2 * c1_exu_coef_MDL1;
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
 * Test Case 3: Simple Solutions
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// This test is made for the Godunov method. The u and w are fixed.

/* ----- ----- ----- ----- ----- -----
 * Coefficients
 * ----- ----- ----- ----- ----- ----- */

// Coefficients used in the model
double pBVal_MDL3;
double c1_exT_coef_MDL3, c2_exT_coef_MDL3, c3_exT_coef_MDL3;
double c2_exu_coef_MDL3, c1_exu_coef_MDL3;
double c2_exw_coef_MDL3, c1_exw_coef_MDL3;

/* ----- ----- ----- ----- ----- -----
 * Domain geometry
 * ----- ----- ----- ----- ----- ----- */

// pB function
double pB_fcn_MDL3(double x) {
	return pBVal_MDL3;
}

// Derivative of pB function
double pBxDer_fcn_MDL3(double x) {
	return 0;
}

/* ----- ----- ----- ----- ----- -----
 * Manufactured solutions / IC
 * ----- ----- ----- ----- ----- ----- */

// Exact T function
double exact_T_fcn_MDL3(double x, double p, double t) {
	return sin(c1_exT_coef_MDL3 * x) * sin(c2_exT_coef_MDL3 * p) * cos(c3_exT_coef_MDL3 * t);
}

// The x-derivative of exact T function
double exact_T_xDer_fcn_MDL3(double x, double p, double t) {
	return c1_exT_coef_MDL3 * cos(c1_exT_coef_MDL3 * x) * sin(c2_exT_coef_MDL3 * p)
			* cos(c3_exT_coef_MDL3 * t);
}

// The p-derivative of exact T function
double exact_T_pDer_fcn_MDL3(double x, double p, double t) {
	return c2_exT_coef_MDL3 * sin(c1_exT_coef_MDL3 * x) * cos(c2_exT_coef_MDL3 * p)
			* cos(c3_exT_coef_MDL3 * t);
}

// The t-derivative of exact T function
double exact_T_tDer_fcn_MDL3(double x, double p, double t) {
	return -c3_exT_coef_MDL3 * sin(c1_exT_coef_MDL3 * x) * sin(c2_exT_coef_MDL3 * p)
			* sin(c3_exT_coef_MDL3 * t);
}

// Exact q function: same with the exact T function
double exact_q_fcn_MDL3(double x, double p, double t) {
	return exact_T_fcn_MDL3(x, p, t);
}

// The x-derivative of exact q function
double exact_q_xDer_fcn_MDL3(double x, double p, double t) {
	return exact_T_xDer_fcn_MDL3(x, p, t);
}

// The p-derivative of exact q function
double exact_q_pDer_fcn_MDL3(double x, double p, double t) {
	return exact_T_pDer_fcn_MDL3(x, p, t);
}

// The t-derivative of exact q function
double exact_q_tDer_fcn_MDL3(double x, double p, double t) {
	return exact_T_tDer_fcn_MDL3(x, p, t);
}

// Exact u function
double exact_u_fcn_MDL3(double x, double p, double t) {
	return 7.5 + cos(c1_exu_coef_MDL3 * x) * cos(c2_exu_coef_MDL3 * p) ;
}

// The x-derivative of exact u function
double exact_u_xDer_fcn_MDL3(double x, double p, double t) {
	return -c1_exu_coef_MDL3 * sin(c1_exu_coef_MDL3 * x) * cos(c2_exu_coef_MDL3 * p) ;
}

// Exact w function
double exact_w_fcn_MDL3(double x, double p, double t) {
	return cos(c1_exw_coef_MDL3 * x) * sin(c2_exw_coef_MDL3 * p) ;
}

// The p derivative exact w function
double exact_w_pDer_fcn_MDL3(double x, double p, double t) {
	return c2_exw_coef_MDL3 * cos(c1_exw_coef_MDL3 * x) * cos(c2_exw_coef_MDL3 * p);
}

/* ----- ----- ----- ----- ----- -----
 * Source functions
 * ----- ----- ----- ----- ----- ----- */

double source_T_fcn_MDL3(double T, double q, double u, double x, double p, double t) {
	return exact_T_tDer_fcn_MDL3(x, p, t)
			+ T * (exact_u_xDer_fcn_MDL3(x, p, t) + exact_w_pDer_fcn_MDL3(x, p, t))
			+ exact_u_fcn_MDL3(x, p, t) * exact_T_xDer_fcn_MDL3(x, p, t)
			+ exact_w_fcn_MDL3(x, p, t) * exact_T_pDer_fcn_MDL3(x, p, t);
}

double source_q_fcn_MDL3(double T, double q, double u, double x, double p, double t) {
	return exact_q_tDer_fcn_MDL3(x, p, t)
			+ q * (exact_u_xDer_fcn_MDL3(x, p, t) + exact_w_pDer_fcn_MDL3(x, p, t))
			+ exact_u_fcn_MDL3(x, p, t) * exact_q_xDer_fcn_MDL3(x, p, t)
			+ exact_w_fcn_MDL3(x, p, t) * exact_q_pDer_fcn_MDL3(x, p, t);
}

double source_u_fcn_MDL3(double T, double q, double u, double x, double p, double t) {
	return 0;
}

void setPar_MDL3() {
	// Parameters on the geometry of the domain
	x0 = 0;
	xf = 100;
	pA = 0;
	pBVal_MDL3 = 100;
	pB_fcnPtr = &pB_fcn_MDL3;
	pBxDer_fcnPtr = &pBxDer_fcn_MDL3;

	int numPeriod = 2;
	c1_exT_coef_MDL3 = numPeriod * TWOPI_CONST / (xf - x0);
	c2_exT_coef_MDL3 = numPeriod * TWOPI_CONST / (pBVal_MDL3 - pA);
	c3_exT_coef_MDL3 = 2;

	c1_exu_coef_MDL3 = 2 * M_PI / xf;
	c2_exu_coef_MDL3 = M_PI / pBVal_MDL3;

	c1_exw_coef_MDL3 = 2 * M_PI / xf;
	c2_exw_coef_MDL3 = 2 * M_PI / pBVal_MDL3;
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Set All Model Parameters
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Select model and set model parameters
void selectModels() {
	// Select model according to modelNo
	if (modelNo == 1) {
		setPar_MDL1();
		return;
	}
	if (modelNo == 2) {
		setPar_MDL2();
		return;
	}
	if (modelNo == 3) {
		setPar_MDL3();
		return;
	}

	// Throw error message when the model number does correspond to any model
	throw "Error: Model does NOT exist!";
}

// Wrapper function to set all parameters for the selected model
void setModels() {
	try {
		selectModels();
	} catch (const char* msg) {
		cerr << msg << endl;
		return;
	}
}

#endif /* MODELS_H_ */
