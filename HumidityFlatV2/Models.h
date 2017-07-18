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

// The velocity u function
double (*uFcnPtr)(double x, double p);
// The velocity omega function
double (*omegaFcnPtr)(double x, double p);

/*
 * The velocity u function
 */
double u_fcn(double x, double p) {
	return (*uFcnPtr)(x, p);
}

/*
 * The velocity omega function
 */
double omega_fcn(double x, double p) {
	return (*omegaFcnPtr)(x, p);
}

/*
 * Cache values for the u and omega functions
 */

double uFcn_cache[numCellsX][numCellsP];
double omegaFcn_cache[numCellsX][numCellsP];

/*
 * Fill in the cache for u and omega, valued at cell centers
 */
void uomega_fillCache() {
	for (int i = 0; i < numCellsX; i++)
		for (int j = 0; j < numCellsP; j++) {
			double x = getCellCenterX(i, j), p = getCellCenterP(i, j);
			uFcn_cache[i][j] = u_fcn(x, p);
			omegaFcn_cache[i][j] = omega_fcn(x, p);
		}
}

/*
 * Return the value of the u function in the (i, j) cell
 */
double getuFcnVal(int i, int j) {
	return uFcn_cache[i][j];
}

/*
 * Return the value of the omega function in the (i, j) cell
 */
double getomegaFcnVal(int i, int j) {
	return omegaFcn_cache[i][j];
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Test Models
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

/* ----- ----- ----- ----- -----
 * Zero initial condition
 * ----- ----- ----- ----- ----- */

double zeroInit(double x, double p, double t, int i, int j) {
	return 0;
}

/* ----- ----- ----- ----- -----
 * Original Model
 * ----- ----- ----- ----- ----- */

/*
 * Preset constants from the original model: constants used in the model functions
 */
const double theta_CONST = 1.25; // Parameter used in the van Leer's minmod limiters
const double p0_CONST = 200;
const double R_CONST = 287;
const double Cp_CONST = 1004;
const double Rv_CONST = 461.50;

/*
 * Coefficients for defining the function u and omega (for efficiency)
 */
const double c1_uomega_fcn_orig = M_PI / p0_CONST,
		c2_uomega_fcn_orig = 2 * M_PI / xL;
/*
 * The function u from the original model
 */
double u_fcn_orig(double x, double p) {
	return 7.5 + cos(c1_uomega_fcn_orig * p) * cos(c2_uomega_fcn_orig * x);
}

/*
 * The function omega from the original model
 */
double omega_fcn_orig(double x, double p) {
	return sin(c1_uomega_fcn_orig * p) * cos(c2_uomega_fcn_orig * x);
}

/*
 * Cache for the derivatives of the u and omega functions from the original model
 */
double uxDer_Orig_cache[numCellsX][numCellsP];
double omegapDer_Orig_cache[numCellsX][numCellsP];

/*
 * Calculate the derivatives of the u and omega functions from the original model
 * and store them in the cache
 */
void uomegaDer_Orig_fillCache() {
	for (int i = 0; i < numCellsX; i++)
		for (int j = 0; j < numCellsP; j++) {
			double x = getCellCenterX(i, j), p = getCellCenterP(i, j);
			double input1 = c1_uomega_fcn_orig * p, input2 = c2_uomega_fcn_orig * x;
			double pTerm_uFcn = cos(input1), xTerm_omegaFcn = cos(input2);
			uxDer_Orig_cache[i][j] = - pTerm_uFcn * c2_uomega_fcn_orig * sin(input2);
			omegapDer_Orig_cache[i][j] = c1_uomega_fcn_orig * cos(input1) * xTerm_omegaFcn;
		}
}

/*
 * Set up the original model
 */
void prep_Orig() {
	uFcnPtr = &u_fcn_orig;
	omegaFcnPtr = &omega_fcn_orig;
	uomega_fillCache();
	uomegaDer_Orig_fillCache();
}

double initT_Orig(double x, double p, double t, int j, int k) {
	return 0;
}

double initq_Orig(double x, double p, double t, int j, int k) {
	return 0;
}

void addSource_Orig(double ans[2], double T, double q, double x, double p, double t, int j, int k) {
	double omegaVal = omegaFcn_cache[j][k];
	double qsVal = 0.622 / p * 6.112 * exp(17.67 * (T - 273.15) / (T - 29.65));
	double deltaVal = 0.25 * (1 - sign(omegaVal)) * (1 + sign(q - qsVal));
	double LVal = 2.5008e6 - 2300 * (T - 275);
	double FVal = qsVal * T * (LVal * R_CONST - Cp_CONST * Rv_CONST * T) / (
			Cp_CONST * Rv_CONST * T * T + qsVal * LVal * LVal);
	ans[0] += omegaVal * (R_CONST * T - deltaVal * LVal * FVal) / (p * Cp_CONST);
	ans[1] += deltaVal * omegaVal * FVal / p;
}

/* ----- ----- ----- ----- -----
 * Test 1
 * ----- ----- ----- ----- ----- */

// Coefficients to be used in the test
const double oscillationFactor_Test1 = 1;
double c1_Test1 = oscillationFactor_Test1 * 4 * M_PI / xL;
double c2_Test1 = oscillationFactor_Test1 * M_PI / 200;
double c3_Test1 = 2 * M_PI;

/*
 * Filling the cache for the test
 */
void prep_Test1() {
	uFcnPtr = &u_fcn_orig;
	omegaFcnPtr = &omega_fcn_orig;
	uomega_fillCache();
	uomegaDer_Orig_fillCache();
}

// Scale variable
double scale_Test1 = 1;

/*
 * Manufactured solution for T
 */
double exact_T_Test1(double x, double p, double t, int j, int k) {
	return sin(c1_Test1 * x) * sin(c2_Test1 * p) * cos(c3_Test1 * t) * scale_Test1;
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
			T * (uxDer_Orig_cache[j][k] + omegapDer_Orig_cache[j][k]) +
			T_tPart * (uFcn_cache[j][k] * c1_Test1 * cos(inputx) * T_pPart +
					omegaFcn_cache[j][k] * c2_Test1 * T_xPart * cos(inputp));
	ans[0] += res1 * scale_Test1;
	ans[1] += 0 * scale_Test1;
}

#endif /* MODELS_H_ */
