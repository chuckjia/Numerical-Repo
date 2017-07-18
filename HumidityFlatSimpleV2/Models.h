/*
 * Models.h
 *
 *  Created on: Jul 14, 2017
 *      Author: chuckjia
 */

#ifndef MODELS_H_
#define MODELS_H_
#include "Mesh.h"

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Velocity Functions (Needed In Manufactured Source Functions)
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

/*
 * Preset constants from the original model: constants used in the model functions
 */
const double p0_CONST = 1000;

/*
 * Coefficients for defining the function u and omega (for efficiency)
 */
const double c1_uomega_fcn_orig = M_PI / p0_CONST,
		c2_uomega_fcn_orig = 2 * M_PI / xL;

/*
 * The velocity u function
 */
double u_fcn(double x, double p) {
	return 7.5 + cos(c1_uomega_fcn_orig * p) * cos(c2_uomega_fcn_orig * x);
}

double u_xDer_fcn(double x, double p) {
	return -c2_uomega_fcn_orig * cos(c1_uomega_fcn_orig * p) * sin(c2_uomega_fcn_orig * x);
}

/*
 * The velocity omega function
 */
double omega_fcn(double x, double p) {
	return sin(c1_uomega_fcn_orig * p) * cos(c2_uomega_fcn_orig * x);
}

double omega_pDer_fcn(double x, double p) {
	return cos(c1_uomega_fcn_orig * p) * cos(c2_uomega_fcn_orig * x);
}

/* ----- ----- ----- ----- -----
 * Test 1
 * ----- ----- ----- ----- ----- */

// Coefficients to be used in the test
double c1_Test1 = 4 * M_PI / 10;
double c2_Test1 = M_PI / 10;
double c3_Test1 = 2 * M_PI;

/*
 * Manufactured solution for T
 */
double exact_T_Test1(double x, double p, double t, int i, int j) {
	return sin(c1_Test1 * x) * sin(c2_Test1 * p) * cos(c3_Test1 * t);
}

/*
 * Manufactured solution for q
 */
double exact_q_Test1(double x, double p, double t, int j, int k) {
	return 0;
}

/*
 * Manufactured source function
 */
void addSource_Test1(double ans[2], double T, double q, double x, double p, double t, int j, int k) {
	double inputx = c1_Test1 * x, inputp = c2_Test1 * p, inputt = c3_Test1 * t;
	double T_xPart = sin(inputx), T_pPart = sin(inputp), T_tPart = cos(inputt);
	double res1 = - c3_Test1 * sin(inputt) * T_xPart * T_pPart +
			T * (u_xDer_fcn(x, p) + omega_pDer_fcn(x, p)) +
			T_tPart * (u_fcn(x, p) * c1_Test1 * cos(inputx) * T_pPart +
					omega_fcn(x,p) * c2_Test1 * T_xPart * cos(inputp));
	ans[0] += res1;
	ans[1] += 0;
}

#endif /* MODELS_H_ */
