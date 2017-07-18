/*
 * Models.h
 *
 *  Created on: Jun 20, 2017
 *      Author: chuckjia
 */

#ifndef MODELS_H_
#define MODELS_H_
#include "Mesh.h"

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Velocity Functions (Needed In Manufactured Sources)
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

double (*uFcnPtr)(double x, double p);
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

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Original Model
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

/*
 * Preset constants from the original model: constants used in the model functions
 */
const double theta_CONST = 1.25; // Parameter used in the van Leer's minmod limiters
const double p0_CONST = 1000;
const double R_CONST = 287;
const double Cp_CONST = 1004;
const double Rv_CONST = 461.50;

/*
 * Coefficients for defining the function u and omega (for efficiency)
 */
const double c1_uomega_fcn_origModel = M_PI / p0_CONST,
		c2_uomega_fcn_origModel = 2 * M_PI / xL;
/*
 * The function u
 */
double u_fcn_orig(double x, double p) {
	return 7.5 + cos(c1_uomega_fcn_origModel * p) * cos(c2_uomega_fcn_origModel * x);
}

/*
 * The function omega
 */
double omega_fcn_orig(double x, double p) {
	return sin(c1_uomega_fcn_origModel * p) * cos(c2_uomega_fcn_origModel * x);
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Test Models
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

/* ----- ----- ----- ----- -----
 * Cache: for efficiency
 * ----- ----- ----- ----- ----- */

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

double getuFcnVal(int i, int j) {
	return uFcn_cache[i][j];
}

double getomegaFcnVal(int i, int j) {
	return omegaFcn_cache[i][j];
}

double uxDer_Orig_cache[numCellsX][numCellsP];
double omegapDer_Orig_cache[numCellsX][numCellsP];

void uomegaDer_Orig_fillCache() {
	for (int i = 0; i < numCellsX; i++)
		for (int j = 0; j < numCellsP; j++) {
			double x = getCellCenterX(i, j), p = getCellCenterP(i, j);
			double input1 = c1_uomega_fcn_origModel * p, input2 = c2_uomega_fcn_origModel * x;
			double pTerm_uFcn = cos(input1), xTerm_omegaFcn = cos(input2);
			uxDer_Orig_cache[i][j] = - pTerm_uFcn * c2_uomega_fcn_origModel * sin(input2);
			omegapDer_Orig_cache[i][j] = c1_uomega_fcn_origModel * cos(input1) * xTerm_omegaFcn;
		}
}

/* ----- ----- ----- ----- -----
 * Zero initial condition
 * ----- ----- ----- ----- ----- */

double zeroInit(double x, double p, double t, int j, int k) {
	return 0;
}

/* ----- ----- ----- ----- -----
 * Test 1
 * ----- ----- ----- ----- ----- */

// Coefficients to be used in the test
const double oscillationFactor = 1;
double c1_Test1 = oscillationFactor * 4 * M_PI / xL;
double c2_Test1 = oscillationFactor * M_PI / 200;
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
double soln_T_Test1(double x, double p, double t, int j, int k) {
	return sin(c1_Test1 * x) * sin(c2_Test1 * p) * cos(c3_Test1 * t) * scale_Test1;
}

/*
 * Manufactured solution for q
 */
double soln_q_Test1(double x, double p, double t, int j, int k) {
	return 0;
}

/*
 * Manufactured source function
 */
void source_Test1(double ans[2], double T, double q, double x, double p, double t, int j, int k) {
	double inputx = c1_Test1 * x, inputp = c2_Test1 * p, inputt = c3_Test1 * t;
	double T_xPart = sin(inputx), T_pPart = sin(inputp), T_tPart = cos(inputt);
	double res1 = - c3_Test1 * sin(inputt) * T_xPart * T_pPart +
			T * (uxDer_Orig_cache[j][k] + omegapDer_Orig_cache[j][k]) +
			T_tPart * (uFcn_cache[j][k] * c1_Test1 * cos(inputx) * T_pPart +
					omegaFcn_cache[j][k] * c2_Test1 * T_xPart * cos(inputp));
	ans[0] += res1 * scale_Test1;
	ans[1] += 0 * scale_Test1;
}

/* ----- ----- ----- ----- -----
 * Test 2
 * ----- ----- ----- ----- ----- */

double sl_cache_Test2[numCellsX][numCellsP][2];

void solnParts_fillCache_Test2() {
	for (int j = 0; j < numCellsX; j++)
		for (int k = 0; k < numCellsP; k++) {
			double x = getCellCenterX(j, k), p = getCellCenterP(j, k);
			sl_cache_Test2[j][k][0] = (x - x0) * (x - xL);
			sl_cache_Test2[j][k][1] = log(p) * (p - pA) * (p - pB);
		}
}

void prep_Test2() {
	uFcnPtr = &u_fcn_orig;  // Need to execute this before prep_Tests
	omegaFcnPtr = &omega_fcn_orig;
	uomega_fillCache();
	solnParts_fillCache_Test2();
}

double coeff_Test2 = 2 * M_PI;
double x0xLSum_Test2 = x0 + xL;
double pApBSum_Test2 = pA + pB;
double scaleCoeff_Test2 = 1e10;

double soln_T_Test2(double x, double p, double t, int j, int k) {
	return cos(coeff_Test2 * t) * sl_cache_Test2[j][k][0] * sl_cache_Test2[j][k][1] / scaleCoeff_Test2;
}

void source_Test2(double ans[2], double T, double q, double x, double p, double t, int j, int k) {
	double solnPartX = sl_cache_Test2[j][k][0], solnPartP = sl_cache_Test2[j][k][1],
			inputt = coeff_Test2 * t, solnPartT = cos(inputt);
	ans[0] += (-coeff_Test2 * sin(inputt) * solnPartX * solnPartP +
			T * (uxDer_Orig_cache[j][k] + omegapDer_Orig_cache[j][k]) +
			uFcn_cache[j][k] * solnPartT * solnPartP * (2 * x - x0xLSum_Test2) +
			omegaFcn_cache[j][k] * solnPartT * solnPartX * (
					(p - pA) * (p - pB) / p + log(p) * (2 * p - pApBSum_Test2))) /
					scaleCoeff_Test2;
	ans[1] += 0;
}

/* ----- ----- ----- ----- -----
 * Test 3: Test 1 with constant u and omega
 * ----- ----- ----- ----- ----- */

// Coefficients to be used in the test
double c1_Test3 = 4 * M_PI / xL;
double c2_Test3 = M_PI / 200;
double c3_Test3 = 2 * M_PI;
double uFcnVal_Test3 = 0;
double omegaFcnVal_Test3 = 0;

double u_fcn_Test3() {
	return uFcnVal_Test3;
}

double omega_fcn_Test3() {
	return omegaFcnVal_Test3;
}

void prep_Test3() {
	uFcnPtr = &u_fcn_Test3;
	omegaFcnPtr = &omega_fcn_Test3;
}

// Scale variable
double scale_Test3 = 1;

/*
 * Manufactured solution for T
 */
double soln_T_Test3(double x, double p, double t, int j, int k) {
	return sin(c1_Test1 * x) * sin(c2_Test1 * p) * cos(c3_Test1 * t) * scale_Test3;
}

/*
 * Manufactured solution for q
 */
double soln_q_Test3(double x, double p, double t, int j, int k) {
	return 0;
}

/*
 * Manufactured source function
 */
void source_Test3(double ans[2], double T, double q, double x, double p, double t, int j, int k) {
	double inputx = c1_Test3 * x, inputp = c2_Test3 * p, inputt = c3_Test3 * t;
	double T_xPart = sin(inputx), T_pPart = sin(inputp), T_tPart = cos(inputt);
	ans[0] = - c3_Test3 * sin(inputt) * T_xPart * T_pPart +
			T_tPart * (uFcnVal_Test3 * c1_Test3 * cos(inputx) * T_pPart +
					omegaFcnVal_Test3 * c2_Test3 * T_xPart * cos(inputp));
	ans[0] += ans[0] * scale_Test3;
	ans[1] += 0;
}

/* ----- ----- ----- ----- -----
 * Test 4
 * ----- ----- ----- ----- ----- */

// Coefficients to be used in the test
double c1_Test4 = 4 * M_PI / (xL - x0);
double c2_Test4 = M_PI / (pB - pA);
double c3_Test4 = 2 * M_PI;

/*
 * Filling the cache for the test
 */
void prep_Test4() {
	uFcnPtr = &u_fcn_orig;
	omegaFcnPtr = &omega_fcn_orig;
	uomega_fillCache();
	uomegaDer_Orig_fillCache();
}


// Scale variable
double scale_Test4 = 1;

/*
 * Manufactured solution for T
 */
double soln_T_Test4(double x, double p, double t, int j, int k) {
	return cos(c1_Test4 * x) * cos(c2_Test4 * p) * cos(c3_Test4 * t) * scale_Test4;
}

/*
 * Manufactured solution for q
 */
double soln_q_Test4(double x, double p, double t, int j, int k) {
	return 0;
}

/*
 * Manufactured source function
 */
void source_Test4(double ans[2], double T, double q, double x, double p, double t, int j, int k) {
	double inputx = c1_Test4 * x, inputp = c2_Test4 * p, inputt = c3_Test4 * t;
	double T_xPart = sin(inputx), T_pPart = sin(inputp), T_tPart = cos(inputt);
	double res1 = - c3_Test4 * sin(inputt) * T_xPart * T_pPart +
			T * (uxDer_Orig_cache[j][k] + omegapDer_Orig_cache[j][k]) -
			T_tPart * (uFcn_cache[j][k] * c1_Test4 * sin(inputx) * T_pPart +
					omegaFcn_cache[j][k] * c2_Test4 * T_xPart * sin(inputp));
	ans[0] += res1 * scale_Test4;
	ans[1] += 0 * scale_Test4;
}

/* ----- ----- ----- ----- -----
 * Test 5
 * ----- ----- ----- ----- ----- */

/*
 * Coefficients for u and omega
 */
double c_uomega = M_PI / (xL - x0);

double u_fcn_Test5(double x, double p) {
	return 7.5 + sin(c_uomega * x);
}

double omega_fcn_Test5(double x, double p) {
	return sin(c_uomega * p);
}

/*
 * Filling the cache for the test
 */
void prep_Test5() {
	uFcnPtr = &u_fcn_Test5;
	omegaFcnPtr = &omega_fcn_Test5;
	uomega_fillCache();
}

// Coefficients to be used in the test
double c1_Test5 = 2 * M_PI / 5;

/*
 * Manufactured solution for T
 */
double soln_T_Test5(double x, double p, double t, int j, int k) {
	return (1 + t * t) * cos(c1_Test5 * x) * cos(c1_Test5 * p);
}

/*
 * Manufactured solution for q
 */
double soln_q_Test5(double x, double p, double t, int j, int k) {
	return 0;
}

/*
 * Manufactured source function
 */
void source_Test5(double ans[2], double T, double q, double x, double p, double t, int j, int k) {
	double inputx = c1_Test5 * x, inputp = c1_Test5 * p;
	double xPart = cos(inputx), pPart = cos(inputp), tPart = 1 + t * t;
	ans[0] += 2 * t * xPart * pPart +
			T * c_uomega * (cos(c_uomega * x) + cos(c_uomega * p)) -
			tPart * c1_Test5 * (sin(c1_Test5 * x) * cos(c1_Test5 * p) * u_fcn_Test5(x, p) +
					cos(c1_Test5 * x) * sin(c1_Test5 * p) * omega_fcn_Test5(x, p));
	ans[1] += 0;
}

/* ----- ----- ----- ----- -----
 * Original Model
 * ----- ----- ----- ----- ----- */

void prep_Orig() {
	uFcnPtr = &u_fcn_orig;  // Need to execute this before prep_Tests
	omegaFcnPtr = &omega_fcn_orig;
}

double initTOrig(double x, double p, double t, int j, int k) {
	return 0;
}

double initqOrig(double x, double p, double t, int j, int k) {
	return 0;
}

void source_Orig(double ans[2], double T, double q, double x, double p, double t, int j, int k) {
	double omegaVal = omegaFcn_cache[j][k];
	double qsVal = 0.622 / p * 6.112 * exp(17.67 * (T - 273.15) / (T - 29.65));
	double deltaVal = 0.25 * (1 - sign(omegaVal)) * (1 + sign(q - qsVal));
	double LVal = 2.5008e6 - 2300 * (T - 275);
	double FVal = qsVal * T * (LVal * R_CONST - Cp_CONST * Rv_CONST * T) / (
			Cp_CONST * Rv_CONST * T * T + qsVal * LVal * LVal);
	ans[0] += omegaVal * (R_CONST * T - deltaVal * LVal * FVal) / (p * Cp_CONST);
	ans[1] += deltaVal * omegaVal * FVal / p;
}

#endif /* MODELS_H_ */
