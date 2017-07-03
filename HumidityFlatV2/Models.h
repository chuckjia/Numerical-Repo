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
 * Original Model
 * ----- ----- ----- ----- ----- */

double initTOrig(double x, double p, double t, int j, int k) {
	return 300 - (1 - p / p0_CONST) * 50;
}

double initqOrig(double x, double p, double t, int j, int k) {
	double T = initTOrig(x, p, t, j, k);
	return 0.622 / p * 6.112 * exp(17.67 * (T - 273.15) / (T - 29.65)) - 0.0052;
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
		// double TVal = exp(t) * T_xpTerms_cache_Test1[j][k];
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
	return cos(coeff3_Test2 * t) * sin(coeff1_Test2 * x) * sin(coeff2_Test2 * p) * 100;
}

void source_Test2(double ans[2], double T, double q, double x, double p, double t, int j, int k) {
	double inputx = coeff1_Test2 * x, inputp = coeff2_Test2 * p, inputt = coeff3_Test2 * t;
	double T_xPart = sin(inputx), T_pPart = sin(inputp), T_tPart = cos(inputt);
	ans[0] = - coeff3_Test2 * sin(inputt) * T_xPart * T_pPart +
			T * (uxDer_cache[j][k] + omegapDer_cache[j][k]) +
			T_tPart * (uFcn_cache[j][k] * coeff1_Test2 * cos(inputx) * T_pPart +
					omegaFcn_cache[j][k] * coeff2_Test2 * T_xPart * cos(inputp));
	ans[0] = ans[0] * 100;
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
double scaleCoeff_Test3 = 1e8;

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

#endif /* MODELS_H_ */
