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
 * Test Models
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

/* ----- ----- ----- ----- -----
 * Cache: for efficiency
 * ----- ----- ----- ----- ----- */

double uFcn_cache[numCellsX][numCellsP];
double omegaFcn_cache[numCellsX][numCellsP];
double uxDer_cache[numCellsX][numCellsP];
double omegapDer_cache[numCellsX][numCellsP];

void u_omega_fillCache() {
	for (int j = 0; j < numCellsX; j++)
		for (int k = 0; k < numCellsP; k++) {
			double x = getCellCenterX(j, k), p = getCellCenterP(j, k);
			double input1 = u_fcn_COEFF1 * p, input2 = omega_fcn_COEFF2 * x;
			double pTerm_uFcn = cos(input1), xTerm_omegaFcn = cos(input2);
			uFcn_cache[j][k] = 7.5 + pTerm_uFcn * cos(input2);
			omegaFcn_cache[j][k] = sin(input1) * xTerm_omegaFcn;
			uxDer_cache[j][k] = - pTerm_uFcn * u_fcn_COEFF2 * sin(input2);
			omegapDer_cache[j][k] = omega_fcn_COEFF1 * cos(input1) * xTerm_omegaFcn;
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

const double x1_Test1 = x0 + (xL - x0) * 0.1;
const double x2_Test1 = x0 + (xL - x0) * 0.9;
const double p1_Test1 = pA + (pB - pA) * 0.1;
const double p2_Test1 = pA + (pB - pA) * 0.9;
const double x1x2Sum_Test1 = x1_Test1 + x2_Test1;
// const double x1x2DiffInv_Test1 = 1 / (x2_Test1 - x1_Test1);
const double x1x2Diff_Test1 = x2_Test1 - x1_Test1;
const double p1p2Sum_Test1 = p1_Test1 + p2_Test1;
// const double p1p2DiffInv_Test1 = 1 / (p2_Test1 - p1_Test1);
const double p1p2Diff_Test1 = p2_Test1 - p1_Test1;

double T_xpTerms_cache_Test1[numCellsX][numCellsP];
double xTerm_cache_Test1[numCellsX][numCellsP];
double pTerm_cache_Test1[numCellsX][numCellsP];

void fillCache_Test1() {
	for (int j = 0; j < numCellsX; j++)
		for (int k = 0; k < numCellsP; k++) {
			double x = getCellCenterX(j, k), p = getCellCenterP(j, k);
			if (x > x1_Test1 && x < x2_Test1 && p > p1_Test1 && p < p2_Test1) {
				double xTerm1 = (2 * x - x1x2Sum_Test1) / x1x2Diff_Test1;
				double pTerm1 = (2 * p - p1p2Sum_Test1) / p1p2Diff_Test1;
				double xTerm = -1 / (1 - xTerm1 * xTerm1);
				double pTerm = -1 / (1 - pTerm1 * pTerm1);
				xTerm_cache_Test1[j][k] = xTerm;
				pTerm_cache_Test1[j][k] = pTerm;
				T_xpTerms_cache_Test1[j][k] = exp(xTerm) * exp(pTerm);
			}
		}
}

void prep_Test1() {
	u_omega_fillCache();
	fillCache_Test1();
}

double soln_T_Test1(double x, double p, double t, int j, int k) {
	if (x < x1_Test1 || x > x2_Test1 || p < p1_Test1 || p > p2_Test1)
		return 0;
	return exp(t) * T_xpTerms_cache_Test1[j][k];
}

void source_Test1(double ans[2], double T, double q, double x, double p, double t, int j, int k) {
	ans[1] = 0;
	if (x < x1_Test1 || x > x2_Test1 || p < p1_Test1 || p > p2_Test1)
		ans[0] = 0;
	else {
		double xTerm = xTerm_cache_Test1[j][k], pTerm = pTerm_cache_Test1[j][k];
		double TxVal_partial = - xTerm * xTerm * 4 * (2 * x - x1x2Sum_Test1) /
				(x1x2Diff_Test1 * x1x2Diff_Test1);
		double TpVal_partial = - pTerm * pTerm * 4 * (2 * p - p1p2Sum_Test1) /
				(p1p2Diff_Test1 * p1p2Diff_Test1);
		ans[0] = T * (1 + uxDer_cache[j][k] + uFcn_cache[j][k] * TxVal_partial +
				omegapDer_cache[j][k] + omegaFcn_cache[j][k] * TpVal_partial);
	}
}

/* ----- ----- ----- ----- -----
 * Test 2
 * ----- ----- ----- ----- ----- */

double coeff1_Test2 = 4 * M_PI / xL;
double coeff2_Test2 = M_PI / 200;
double coeff3_Test2 = 2 * M_PI;

void prep_Test2() {
	u_omega_fillCache();
}

double soln_T_Test2(double x, double p, double t, int j, int k) {
	return cos(coeff3_Test2 * t) * sin(coeff1_Test2 * x) * sin(coeff2_Test2 * p) * 1000;
}

void source_Test2(double ans[2], double T, double q, double x, double p, double t, int j, int k) {
	double inputx = coeff1_Test2 * x, inputp = coeff2_Test2 * p, inputt = coeff3_Test2 * t;
	double T_xPart = sin(inputx), T_pPart = sin(inputp), T_tPart = cos(inputt);
	ans[0] = - coeff3_Test2 * sin(inputt) * T_xPart * T_pPart +
			T * (uxDer_cache[j][k] + omegapDer_cache[j][k]) +
			T_tPart * (uFcn_cache[j][k] * coeff1_Test2 * cos(inputx) * T_pPart +
					omegaFcn_cache[j][k] * coeff2_Test2 * T_xPart * cos(inputp));
	ans[0] = ans[0] * 1000;
	ans[1] = 0;
}

/* ----- ----- ----- ----- -----
 * Test 3
 * ----- ----- ----- ----- ----- */

double solnCache_Test3[numCellsX][numCellsP][2];

void solnParts_fillCache_Test3() {
	for (int j = 0; j < numCellsX; j++)
		for (int k = 0; k < numCellsP; k++) {
			double x = getCellCenterX(j, k), p = getCellCenterP(j, k);
			solnCache_Test3[j][k][0] = (x - x0) * (x - xL);
			solnCache_Test3[j][k][1] = log(p) * (p - pA) * (p - pB);
		}
}

void prep_Test3() {
	u_omega_fillCache();
	solnParts_fillCache_Test3();
}

double coeff_Test3 = 2 * M_PI;
double x0xLSum_Test3 = x0 + xL;
double pApBSum_Test3 = pA + pB;
double scaleCoeff_Test3 = 1e10;

double soln_T_Test3(double x, double p, double t, int j, int k) {
	return cos(coeff_Test3 * t) * solnCache_Test3[j][k][0] * solnCache_Test3[j][k][1] / scaleCoeff_Test3;
}

void source_Test3(double ans[2], double T, double q, double x, double p, double t, int j, int k) {
	double solnPartX = solnCache_Test3[j][k][0], solnPartP = solnCache_Test3[j][k][1],
			inputt = coeff_Test3 * t, solnPartT = cos(inputt);
	ans[0] = (-coeff_Test3 * sin(inputt) * solnPartX * solnPartP +
			T * (uxDer_cache[j][k] + omegapDer_cache[j][k]) +
			uFcn_cache[j][k] * solnPartT * solnPartP * (2 * x - x0xLSum_Test3) +
			omegaFcn_cache[j][k] * solnPartT * solnPartX * (
					(p - pA) * (p - pB) / p + log(p) * (2 * p - pApBSum_Test3))) /
					scaleCoeff_Test3;
	ans[1] = 0;
}

/* ----- ----- ----- ----- -----
 * Test 4
 * ----- ----- ----- ----- ----- */

void prep_Test4() {
}

double soln_T_Test4(double x, double p, double t, int j, int k) {
	return 1;
}

void source_Test4(double ans[2], double T, double q, double x, double p, double t, int j, int k) {
	ans[0] = 0;
	ans[1] = 0;
}

/* ----- ----- ----- ----- -----
 * Test 5: Test 2 with omegaFcn == 0
 * ----- ----- ----- ----- ----- */

double coeff1_Test5 = 4 * M_PI / xL;
double coeff2_Test5 = M_PI / 200;
double coeff3_Test5 = 2 * M_PI;

void prep_Test5() {
	u_omega_fillCache();
}

double soln_T_Test5(double x, double p, double t, int j, int k) {
	return cos(coeff3_Test5 * t) * sin(coeff1_Test5 * x) * sin(coeff2_Test5 * p) * 1000;
}

void source_Test5(double ans[2], double T, double q, double x, double p, double t, int j, int k) {
	double inputx = coeff1_Test5 * x, inputp = coeff2_Test5 * p, inputt = coeff3_Test5 * t;
	double T_xPart = sin(inputx), T_pPart = sin(inputp), T_tPart = cos(inputt);
	ans[0] = - coeff3_Test5 * sin(inputt) * T_xPart * T_pPart +
			T * uxDer_cache[j][k] +
			T_tPart * (uFcn_cache[j][k] * coeff1_Test5 * cos(inputx) * T_pPart);
	ans[0] = ans[0] * 1000;
	ans[1] = 0;
}

/* ----- ----- ----- ----- -----
 * Test 6: Test 2 with uFcn == 0
 * ----- ----- ----- ----- ----- */

double coeff1_Test6 = 4 * M_PI / xL;
double coeff2_Test6 = M_PI / 200;
double coeff3_Test6 = 2 * M_PI;

void prep_Test6() {
	u_omega_fillCache();
}

double soln_T_Test6(double x, double p, double t, int j, int k) {
	return cos(coeff3_Test6 * t) * sin(coeff1_Test6 * x) * sin(coeff2_Test6 * p) * 1000;
}

void source_Test6(double ans[2], double T, double q, double x, double p, double t, int j, int k) {
	double inputx = coeff1_Test6 * x, inputp = coeff2_Test6 * p, inputt = coeff3_Test6 * t;
	double T_xPart = sin(inputx), T_pPart = sin(inputp), T_tPart = cos(inputt);
	ans[0] = - coeff3_Test6 * sin(inputt) * T_xPart * T_pPart +
			T * omegapDer_cache[j][k] +
			T_tPart * (omegaFcn_cache[j][k] * coeff1_Test6 * cos(inputx) * T_pPart);
	ans[0] = ans[0] * 1000;
	ans[1] = 0;
}

/* ----- ----- ----- ----- -----
 * Test 7: Test 2 with uFcn == 1 and omegaFcn == 1
 * ----- ----- ----- ----- ----- */

double coeff1_Test7 = 4 * M_PI / xL;
double coeff2_Test7 = M_PI / 200;
double coeff3_Test7 = 2 * M_PI;

void prep_Test7() {
}

double soln_T_Test7(double x, double p, double t, int j, int k) {
	return cos(coeff3_Test7 * t) * sin(coeff1_Test7 * x) * sin(coeff2_Test7 * p);
}

void source_Test7(double ans[2], double T, double q, double x, double p, double t, int j, int k) {
	double inputx = coeff1_Test7 * x, inputp = coeff2_Test7 * p, inputt = coeff3_Test7 * t;
	double T_xPart = sin(inputx), T_pPart = sin(inputp), T_tPart = cos(inputt);
	ans[0] = - coeff3_Test7 * sin(inputt) * T_xPart * T_pPart +
			T_tPart * (T_xPart + T_pPart);
	ans[1] = 0;
}

/* ----- ----- ----- ----- -----
 * Test 8: Test 2 without cache values
 * ----- ----- ----- ----- ----- */

double coeff1_Test8 = 4 * M_PI / xL;
double coeff2_Test8 = M_PI / 200;
double coeff3_Test8 = 2 * M_PI;

void prep_Test8() {
}

double soln_T_Test8(double x, double p, double t, int j, int k) {
	return cos(coeff3_Test8 * t) * sin(coeff1_Test8 * x) * sin(coeff2_Test8 * p) * 10;
}

void source_Test8(double ans[2], double T, double q, double x, double p, double t, int j, int k) {
	double uxDerVal = - u_fcn_COEFF2 * cos(u_fcn_COEFF1 * p) * sin(u_fcn_COEFF2 * x),
			omegapDerVal = u_fcn_COEFF1 * cos(u_fcn_COEFF1 * p) * cos(u_fcn_COEFF2 * x);
	double TxVal = coeff1_Test8 * cos(coeff3_Test8 * t) * cos(coeff1_Test8 * x) * sin(coeff2_Test8 * p),
			TpVal = coeff2_Test8 * cos(coeff3_Test8 * t) * sin(coeff1_Test8 * x) * cos(coeff2_Test8 * p);
	ans[0] = -coeff3_Test8 * sin(coeff3_Test8 * t) * sin(coeff1_Test8 * x) *
			sin(coeff2_Test8 * p) + T * (uxDerVal + omegapDerVal) +
			u_fcn(x, p) * TxVal + omega_fcn(x, p) * TpVal;
	ans[0] = ans[0] * 10;
	ans[1] = 0;
}

/* ----- ----- ----- ----- -----
 * Test 9: Test 2 with u == 0 and omega == 0
 * ----- ----- ----- ----- ----- */

double coeff1_Test9 = 4 * M_PI / xL;
double coeff2_Test9 = M_PI / 200;
double coeff3_Test9 = 2 * M_PI;

void prep_Test9() {
}

double soln_T_Test9(double x, double p, double t, int j, int k) {
	return cos(coeff3_Test9 * t) * sin(coeff1_Test9 * x) * sin(coeff2_Test9 * p) * 1000;
}

void source_Test9(double ans[2], double T, double q, double x, double p, double t, int j, int k) {
	ans[0] = -coeff3_Test9 * sin(coeff3_Test9 * t) * sin(coeff1_Test9 * x) *
			sin(coeff2_Test9 * p);
	ans[0] = ans[0] * 1000;
	ans[1] = 0;
}

/* ----- ----- ----- ----- -----
 * Test 10: Test 2 with u == const and omega == const
 * ----- ----- ----- ----- ----- */

double coeff1_Test10 = 4 * M_PI / xL;
double coeff2_Test10 = M_PI / 200;
double coeff3_Test10 = 2 * M_PI;

void prep_Test10() {
}

double soln_T_Test10(double x, double p, double t, int j, int k) {
	return cos(coeff3_Test10 * t) * sin(coeff1_Test10 * x) * sin(coeff2_Test10 * p) * 1000;
}

void source_Test10(double ans[2], double T, double q, double x, double p, double t, int j, int k) {
	ans[0] = -coeff3_Test10 * sin(coeff3_Test10 * t) * sin(coeff1_Test10 * x) *
			sin(coeff2_Test10 * p) +
			u_fcn(x, p) * cos(coeff3_Test10 * t) * coeff1_Test10 *
			cos(coeff1_Test10 * x) * sin(coeff2_Test10 * p) +
			omega_fcn(x, p) * cos(coeff3_Test10 * t) * sin(coeff1_Test10 * x) *
			coeff2_Test10 * cos(coeff2_Test10 * p);
	ans[0] = ans[0] * 1000;
	ans[1] = 0;
}

/* ----- ----- ----- ----- -----
 * Original Model
 * ----- ----- ----- ----- ----- */

void prep_Orig() {
	fillCache_Test1();
}

double initTOrig(double x, double p, double t, int j, int k) {
	return soln_T_Test1(x, p, t, j, k);
}

double initqOrig(double x, double p, double t, int j, int k) {
	return 0;
	// return initTOrig(x, p, t, j, k);
}

void source_Orig(double ans[2], double T, double q, double x, double p, double t, int j, int k) {
	double omegaVal = omega_fcn(x, p);
	double qsVal = 0.622 / p * 6.112 * exp(17.67 * (T - 273.15) / (T - 29.65));
	double deltaVal = 0.25 * (1 - sign(omegaVal)) * (1 + sign(q - qsVal));
	double LVal = 2.5008e6 - 2300 * (T - 275);
	double FVal = qsVal * T * (LVal * R_CONST - Cp_CONST * Rv_CONST * T) / (
			Cp_CONST * Rv_CONST * T * T + qsVal * LVal * LVal);
	ans[0] = omegaVal * (R_CONST * T - deltaVal * LVal * FVal) / (p * Cp_CONST);
	ans[1] = deltaVal * omegaVal * FVal / p;
}

#endif /* MODELS_H_ */
