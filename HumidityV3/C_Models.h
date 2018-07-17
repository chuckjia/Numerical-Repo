/*
 * C_Models.h
 *
 *  Created on: Oct 14, 2017
 *      Author: chuckjia
 *
 *  This file contains functions and global variables that determines different PDE models.
 */

#ifndef C_MODELS_H_
#define C_MODELS_H_
#include "B_Constants.h"

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Geometry of the Domain
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

double x0, xf, pA;
double (*pB_fptr)(double x);  // Function pointer: pB function
double (*pBxDer_fptr)(double x);  // Function pointer: derivative of pB function


/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Model 0: Physical Model
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

/* ----- ----- ----- ----- ----- -----
 * Coefficients
 * ----- ----- ----- ----- ----- ----- */

// Coefficients used in the model
double _h_pB_MDL0, _c1_pBxDer_MDL0, _n_initU_MDL0, _cp_initU_MDL0, _cx_initU_MDL0;

// Set values to the coefficients defined in this section
// Coefficients need to be initialized here, because some of them depend on x0, xf, etc
void setFcnCoef_MDL0() {
	_h_pB_MDL0 = 250.;  // Multiplicative factor that controls the height of mountain
	_c1_pBxDer_MDL0 = _h_pB_MDL0 / 18e6;  // heightFactor * 2 / 6000^2

	_n_initU_MDL0 = 1.;   // n in (4.10)
	_cp_initU_MDL0 = M_PI * p0Inv_CONST;   // Coefficient for p in (4.10)
	_cx_initU_MDL0 = 2 * _n_initU_MDL0 * M_PI / xf;  // Coefficient for x in (4.10)
}

/* ----- ----- ----- ----- ----- -----
 * Domain geometry
 * ----- ----- ----- ----- ----- ----- */

// pB function
double pB_fcn_MDL0(double x) {
	double term = x - 37500.;
	return 1000. - _h_pB_MDL0 * exp(-term * term / 36e6);
}

// Derivative of pB function
double pBxDer_fcn_MDL0(double x) {
	double term = x - 37500.;
	return _c1_pBxDer_MDL0 * term * exp(-term * term / 36e6);
}

/* ----- ----- ----- ----- ----- -----
 * Initial Conditions
 * ----- ----- ----- ----- ----- ----- */

// Initial T function
double init_T_fcn_MDL0(double x, double p, double t) {
	return T0_CONST - (1 - p * p0Inv_CONST) * DeltaT_CONST;
}

// Initial q function
double init_q_fcn_MDL0(double x, double p, double t) {
	double T = init_T_fcn_MDL0(x, p, 0);
	return qs_fcn(T, p) - 0.000052;
}

// Initial u function
double init_u_fcn_MDL0(double x, double p, double t) {
	return 7.5 + 2 * cos(_cp_initU_MDL0 * p) * cos(_cx_initU_MDL0 * x);
}

/* ----- ----- ----- ----- ----- -----
 * Source solutions
 * ----- ----- ----- ----- ----- ----- */

// Source function for the T equation
double source_T_fcn_MDL0(double T, double q, double u, double w, double x, double p, double t) {
	double qsVal = qs_fcn(T, p),
			deltaVal = delta_fcn(q, w, qsVal),
			LVal = L_fcn(T),
			FVal = F_fcn(T, qsVal, LVal);
	return w / (p * Cp_CONST) * (R_CONST * T - deltaVal * LVal * FVal);
}

// Source function for the q equation
double source_q_fcn_MDL0(double T, double q, double u, double w, double x, double p, double t) {
	double qsVal = qs_fcn(T, p),
			deltaVal = delta_fcn(q, w, qsVal),
			LVal = L_fcn(T),
			FVal = F_fcn(T, qsVal, LVal);
	return deltaVal * FVal * w / p;
}

// Source function for the u equation
double source_u_fcn_MDL0(double T, double q, double u, double w, double x, double p, double t) {
	return 0;
}

// Set all parameters in Model 1
void setPar_MDL0() {
	// Parameters for the domain geometry
	x0 = 0.;
	xf = 75000.;
	pA = 250.;
	pB_fptr = &pB_fcn_MDL0;
	pBxDer_fptr = &pBxDer_fcn_MDL0;
	// Function coefficients
	setFcnCoef_MDL0();
}


/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Test Case 1: Test Case from the Section 4.1
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

/* ----- ----- ----- ----- ----- -----
 * Coefficients
 * ----- ----- ----- ----- ----- ----- */

// Coefficients used in the model
double _xfCubed_MDL1, _xf6thPow_MDL1, _pBx0_MDL1;
double _c1_exT_MDL1, _c2_exT_MDL1, _c3_exT_MDL1, _c4_exT_MDL1;
double _c1_exU_MDL1, _c1_exUpDer_MDL1, _c1_exW_MDL1;

void setFcnCoef_MDL1() {
	_xfCubed_MDL1 = xf * xf * xf;
	_xf6thPow_MDL1 = pow(xf, 6);
	_pBx0_MDL1 = (*pB_fptr)(x0);  // pB(x0)

	_c1_exT_MDL1 = -1. / (g_CONST * _xfCubed_MDL1);
	_c2_exT_MDL1 = g_CONST / (150 * 450 * 450 * R_CONST);
	_c3_exT_MDL1 = DeltaT_CONST / p0_CONST;
	_c4_exT_MDL1 = DeltaT_CONST - T0_CONST;

	_c1_exU_MDL1 = -3e-12 / _xf6thPow_MDL1;
	_c1_exUpDer_MDL1 = 2 * _c1_exU_MDL1;
	_c1_exW_MDL1 = 3e-12 / _xf6thPow_MDL1;
}

/* ----- ----- ----- ----- ----- -----
 * Domain geometry
 * ----- ----- ----- ----- ----- ----- */

// pB function
double pB_fcn_MDL1(double x) {
	double term = x - 25000.;
	return 1000. - 200. * exp(-term * term / 9e6);
}

// The derivative of pB function
double pBxDer_fcn_MDL1(double x) {
	double term = x - 25000.;
	return term * exp(-term * term / 9e6) / 22500;
}

/* ----- ----- ----- ----- ----- -----
 * Manufactured solutions / IC
 * ----- ----- ----- ----- ----- ----- */

// Exact T function excluding the t-term
double exact_T_helper_fcn_MDL1(double x, double p) {
	return _c1_exT_MDL1 * x * pow(x - xf, 2) * (
			p * (_c2_exT_MDL1 * pow(p - pB_fcn_MDL1(x), 2) - _c3_exT_MDL1) + _c4_exT_MDL1
	);
}

// Manufactured solution: exact T function
double exact_T_fcn_MDL1(double x, double p, double t) {
	return cos(TWO_PI * t) * exact_T_helper_fcn_MDL1(x, p);
}

// x-derivative of the exact T function
double exact_TxDer_fcn_MDL1(double x, double p, double t) {
	double x_minus_xf = x - xf, pB_minus_p = pB_fcn_MDL1(x) - p;
	double xPart1 = x * x_minus_xf * x_minus_xf,
			xPart2 = p * (_c2_exT_MDL1 * pB_minus_p * pB_minus_p - _c3_exT_MDL1) + _c4_exT_MDL1;
	return _c1_exT_MDL1 * cos(TWO_PI * t) * (
			(3 * x - xf) * x_minus_xf * xPart2 +
			xPart1 * (2 * _c2_exT_MDL1 * pB_minus_p * pBxDer_fcn_MDL1(x) * p)
	);
}

// The p-derivative of the exact T function
double exact_TpDer_fcn_MDL1(double x, double p, double t) {
	double p_minus_pB = p - pB_fcn_MDL1(x);
	return _c1_exT_MDL1 * cos(TWO_PI * t) * x * pow(x - xf, 2) * (
			_c2_exT_MDL1 * p_minus_pB * p_minus_pB - _c3_exT_MDL1 +
			p * _c2_exT_MDL1 * 2 * p_minus_pB
	);
}

// The t-derivative of the exact T function
double exact_TtDer_fcn_MDL1(double x, double p, double t) {
	return -TWO_PI * sin(TWO_PI * t) * exact_T_helper_fcn_MDL1(x, p);
}

// Manufactured solution: exact q function
double exact_q_fcn_MDL1(double x, double p, double t) {
	return 0;
}

// The terms in u that involves only x and p
double exact_u_fcn_helper_MDL1(double x, double p) {
	double p_minus_pA = p - pA, p_minus_pB = p - pB_fcn_MDL1(x),
			xTerm = x * (x - xf), pTerm = p_minus_pA * p_minus_pB;
	return _c1_exU_MDL1 * xTerm * xTerm * xTerm * pTerm * pTerm * (p_minus_pA + p_minus_pB);
}

// Manufactured solution: exact u function
double exact_U_fcn_MDL1(double x, double p, double t) {
	return (cos(TWO_PI * t) + 20) * exact_u_fcn_helper_MDL1(x, p);
}

// The x-derivative of the exact u function
double exact_UxDer_fcn_MDL1(double x, double p, double t) {
	double x_minus_xf = x - xf,
			p_minus_pA = p - pA, pB_minus_p = pB_fcn_MDL1(x) - p,
			pB_xDer_val = pBxDer_fcn_MDL1(x);
	double xPart1 = x * x * x, xPart2 = x_minus_xf * x_minus_xf * x_minus_xf,
			xPart3 = pB_minus_p * pB_minus_p, xPart4 = p_minus_pA - pB_minus_p;
	return _c1_exU_MDL1 * p_minus_pA * p_minus_pA * (cos(TWO_PI * t) + 20) * (
			3 * x * x * xPart2 * xPart3 * xPart4
			+ xPart1 * 3 * x_minus_xf * x_minus_xf * xPart3 * xPart4
			+ xPart1 * xPart2 * 2 * pB_minus_p * pB_xDer_val * xPart4
			- xPart1 * xPart2 * xPart3 * pB_xDer_val
	);
}

// The p-derivative of the exact u function
double exact_UpDer_fcn_MDL1(double x, double p, double t) {
	double x_minus_xf = x - xf,
			p_minus_pA = p - pA, p_minus_pB = p - pB_fcn_MDL1(x);
	double p_sum_part = p_minus_pA + p_minus_pB, p_prod_part = p_minus_pA * p_minus_pB;
	return _c1_exUpDer_MDL1 * pow(x * x_minus_xf, 3) * (cos(TWO_PI * t) + 20)
			* p_prod_part * (p_sum_part * p_sum_part + p_prod_part);
}

// The t-derivative of the exact u function
double exact_UtDer_fcn_MDL1(double x, double p, double t) {
	return -TWO_PI * sin(TWO_PI * t) * exact_u_fcn_helper_MDL1(x, p);
}

// Manufactured solution: exact w function
double exact_W_fcn_MDL1(double x, double p, double t) {
	double xPart1 = x * (x - xf), p_minus_pB = p - pB_fcn_MDL1(x);
	return _c1_exW_MDL1 * pow(p - pA, 3) * (cos(TWO_PI * t) + 20)
			* xPart1 * xPart1 * p_minus_pB * p_minus_pB * (
					p_minus_pB * (2 * x - xf) - pBxDer_fcn_MDL1(x) * xPart1
			);
}

// The p-derivative of the exact w function
double exact_WpDer_fcn_MDL1(double x, double p, double t) {
	double xTerm = x * (x - xf), p_minus_pB = p - pB_fcn_MDL1(x), p_minus_pA = p - pA,
			two_x_minus_xf = 2 * x - xf;
	double pPart1 = pow(p_minus_pA, 3), pPart2 = p_minus_pB * p_minus_pB,
			pPart3 = p_minus_pB * two_x_minus_xf - pBxDer_fcn_MDL1(x) * xTerm;
	return _c1_exW_MDL1 * xTerm * xTerm * (cos(TWO_PI * t) + 20) * (
			3 * p_minus_pA * p_minus_pA * pPart2 * pPart3
			+ pPart1 * 2 * p_minus_pB * pPart3
			+ pPart1 * pPart2 * two_x_minus_xf
	);
}

/* ----- ----- ----- ----- ----- -----
 * Source solutions
 * ----- ----- ----- ----- ----- ----- */

// Source function for the T equation
double source_T_fcn_MDL1(double T, double q, double u, double w, double x, double p, double t) {
	return exact_TtDer_fcn_MDL1(x, p, t)
			+ exact_U_fcn_MDL1(x, p, t) * exact_TxDer_fcn_MDL1(x, p, t)
			+ exact_W_fcn_MDL1(x, p, t) * exact_TpDer_fcn_MDL1(x, p, t);
}

// Source function for the q equation
double source_q_fcn_MDL1(double T, double q, double u, double w, double x, double p, double t) {
	return 0;
}

// Source function for the u equation
double source_u_fcn_MDL1(double T, double q, double u, double w, double x, double p, double t) {
	return exact_UtDer_fcn_MDL1(x, p, t)
			+ exact_U_fcn_MDL1(x, p, t) * exact_UxDer_fcn_MDL1(x, p, t)
			+ exact_W_fcn_MDL1(x, p, t) * exact_UpDer_fcn_MDL1(x, p, t);
}

// Set all parameters in model 1
void setPar_MDL1() {
	// Parameters on the domain geometry
	x0 = 0.;
	xf = 50000.;
	pA = 200.;
	pB_fptr = &pB_fcn_MDL1;
	pBxDer_fptr = &pBxDer_fcn_MDL1;
	// Set function coefficients
	setFcnCoef_MDL1();
}


/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Set All
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Wrapper function to set all parameters for the selected model
void setModels() {
	try {
		// Select model according to modelNo
		switch (modelNo) {
		case 0:
			setPar_MDL0();
			return;
		case 1:
			setPar_MDL1();
			return;
		default: // Throw error message when the model number does correspond to any model
			throw "Error: Model does NOT exist!";
		}
	} catch (const char* msg) {
		cerr << msg << endl;
		exit(EXIT_FAILURE);
	}
}

#endif /* C_MODELS_H_ */
