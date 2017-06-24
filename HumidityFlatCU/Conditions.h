/*
 * Conditions.h
 *
 *  Created on: Jun 21, 2017
 *      Author: chuckjia
 */

#ifndef CONDITIONS_H_
#define CONDITIONS_H_
#include "Mesh.h"

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * All Tests Declarations
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Original model
void sourceFcnOrig(double ans[2], double T, double q, double x, double p, double t);
// Test 1
double soln_T_test1(double x, double p, double t);
double source1_test1(double T, double q, double x, double p, double t);

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Initial Conditions
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

/*
 * Initialize the solution
 */
double sl[numCellsX][numCellsP][2];

/*
 * Set initial conditions
 */
void setInitCond() {
	for (int i = 0; i < numCellsX; i++)
		for (int j = 0; j < numCellsP; j++) {
			double x = getCellCenterX(i, j), p = getCellCenterP(i, j);
			sl[i][j][0] = soln_T_test1(x, p, 0);
			sl[i][j][1] = 0;
			// printf("x = %f, p = %f, T = %f\n", x, p, sl[i][j][0]);
		}
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Source Terms
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

/*
 * Source function: wrapper
 */
void calcSourceFcn(double ans[2], double T, double q, double x, double p, double t) {
	ans[0] = source1_test1(T, q, x, p, t);
	ans[1] = 0;
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Test Models
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

/* ----- ----- ----- ----- -----
 * Original Model
 * ----- ----- ----- ----- ----- */

void sourceFcnOrig(double ans[2], double T, double q, double x, double p, double t) {
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

const double x1_Test1 = x0 + (xL - x0) * 0.2;
const double x2_Test1 = x0 + (xL - x0) * 0.8;
const double p1_Test1 = pA + (pB - pA) * 0.2;
const double p2_Test1 = pA + (pB - pA) * 0.8;
const double x1x2Sum_Test1 = x1_Test1 + x2_Test1;
const double x1x2DiffInv_Test1 = 1 / (x2_Test1 - x1_Test1);
const double p1p2Sum_Test1 = p1_Test1 + p2_Test1;
const double p1p2DiffInv_Test1 = 1 / (p2_Test1 - p1_Test1);

double soln_T_test1(double x, double p, double t) {
	double xTerm = (2 * x - x1x2Sum_Test1) * x1x2DiffInv_Test1;
	double pTerm = (2 * p - p1p2Sum_Test1) * p1p2DiffInv_Test1;
	return exp(t) * exp(-1 / (1 - xTerm * xTerm)) * exp(-1 / (1 - pTerm * pTerm));
}

double source1_test1(double T, double q, double x, double p, double t) {
	double tTerm = exp(t);
	double xTerm1 = pow((2 * x - x1x2Sum_Test1) * x1x2DiffInv_Test1, 2),
			xTerm2 = 1 / (1 - xTerm1);
	double pTerm1 = pow((2 * p - p1p2Sum_Test1) * p1p2DiffInv_Test1, 2),
			pTerm2 = 1 / (1 - pTerm1);
	double TVal = tTerm * exp(-xTerm2) * exp(-pTerm2);
	double TxVal = TVal * xTerm2 * xTerm2 * 4 * (2 * x - x1x2Sum_Test1) * x1x2DiffInv_Test1 * x1x2DiffInv_Test1;
	double TpVal = TVal * pTerm2 * pTerm2 * 4 * (2 * p - p1p2Sum_Test1) * p1p2DiffInv_Test1 * p1p2DiffInv_Test1;
	// Define u, u_x, omega, and omega_p
	double input1 = u_fcn_COEFF1 * p, input2 = omega_fcn_COEFF2 * x;
	double pTerm_uFcn =  cos(input1), xTerm_omegaFcn = cos(input2);
	double uVal = 7.5 + pTerm_uFcn * cos(input2),
			omegaVal = sin(input1) * xTerm_omegaFcn;
	double uxVal = - pTerm_uFcn * sin(input2),
			omegapVal = cos(input1) * xTerm_omegaFcn;
	return TVal + uxVal * TVal + uVal * TxVal + omegapVal * TVal + omegaVal * TpVal;
}

#endif /* CONDITIONS_H_ */
