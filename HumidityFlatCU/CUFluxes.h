/*
 * Fluxes.h
 *
 *  Created on: Jun 22, 2017
 *      Author: chuckjia
 */

#ifndef CUFLUXES_H_
#define CUFLUXES_H_
#include "Conditions.h"

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Initialize the Fluxes
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Initialize the fluxes
// Hx[j][k] represents H_{j+1/2,k}^x and Hp[j][k] represents H_{j,k+1/2}^y
double Hx[numCellsX][numCellsP][2];
double Hp[numCellsX][numCellsP][2];

// Initialize parameters for the reconstructed solutions
double sl_x[numCellsX][numCellsP][2];
double sl_p[numCellsX][numCellsP][2];

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Central Upwind: Calculate the Fluxes
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

/*
 * Reconstruct the solution
 */
void rcstrSolnCU() {
	for (int j = 1; j < lastIndexX; j++)
		for (int k = 1; k < lastIndexP; k++) {
			double slThis1 = sl[j][k][0], slThis2 = sl[j][k][1];
			double slLeft1 = sl[j - 1][k][0], slLeft2 = sl[j - 1][k][1],
					slRight1 = sl[j + 1][k][0], slRight2 = sl[j + 1][k][1];
			double slLower1 = sl[j][k - 1][0], slLower2 = sl[j][k - 1][1],
					slUpper1 = sl[j][k + 1][0], slUpper2 = sl[j][k + 1][1];
			sl_x[j][k][0] = minmod3(theta_CONST * (slRight1 - slThis1) * DxInv,
					0.5 * (slRight1 - slLeft1) * DxInv,
					theta_CONST * (slThis1 - slLeft1) * DxInv);
			sl_x[j][k][1] = minmod3(theta_CONST * (slRight2 - slThis2) * DxInv,
					0.5 * (slRight2 - slLeft2) * DxInv,
					theta_CONST * (slThis2 - slLeft2) * DxInv);
			sl_p[j][k][0] = minmod3(theta_CONST * (slUpper1 - slThis1) * DpInv,
					0.5 * (slUpper1 - slLower1) * DpInv,
					theta_CONST * (slThis1 - slLower1) * DpInv);
			sl_p[j][k][1] = minmod3(theta_CONST * (slUpper2 - slThis2) * DpInv,
					0.5 * (slUpper2 - slLower2) * DpInv,
					theta_CONST * (slThis2 - slLower2) * DpInv);
		}

	// When j = 0 or lastIndexX
	int indexRange[2] = {0, lastIndexX};
	for (int ii = 0; ii < 2; ii++) {
		int j = indexRange[ii];
		for (int k = 1; k < lastIndexP; k++) {
			double slThis1 = sl[j][k][0], slThis2 = sl[j][k][1];
			double slLower1 = sl[j][k - 1][0], slLower2 = sl[j][k - 1][1],
					slUpper1 = sl[j][k + 1][0], slUpper2 = sl[j][k + 1][1];
			sl_p[j][k][0] = minmod3(theta_CONST * (slUpper1 - slThis1) * DpInv,
					0.5 * (slUpper1 - slLower1) * DpInv,
					theta_CONST * (slThis1 - slLower1) * DpInv);
			sl_p[j][k][1] = minmod3(theta_CONST * (slUpper2 - slThis2) * DpInv,
					0.5 * (slUpper2 - slLower2) * DpInv,
					theta_CONST * (slThis2 - slLower2) * DpInv);
		}
	}

	// When k = 0 or lastIndexP
	indexRange[1] = lastIndexP;
	for (int ii = 0; ii < 2; ii++) {
		int k = indexRange[ii];
		for (int j = 0; j < lastIndexX; j++) {
			double slThis1 = sl[j][k][0], slThis2 = sl[j][k][1];
			double slLeft1 = sl[j - 1][k][0], slLeft2 = sl[j - 1][k][1],
					slRight1 = sl[j + 1][k][0], slRight2 = sl[j + 1][k][1];
			sl_x[j][k][0] = minmod3(theta_CONST * (slRight1 - slThis1) * DxInv,
					0.5 * (slRight1 - slLeft1) * DxInv,
					theta_CONST * (slThis1 - slLeft1) * DxInv);
			sl_x[j][k][1] = minmod3(theta_CONST * (slRight2 - slThis2) * DxInv,
					0.5 * (slRight2 - slLeft2) * DxInv,
					theta_CONST * (slThis2 - slLeft2) * DxInv);
		}
	}
}

/*
 * Calculate the reconstructed solution value
 */
void calcRcSolnCU(double ans[2], int j, int k, double x, double p) {
	double xCenter = getCellCenterX(j, k), pCenter = getCellCenterP(j, k);
	double xTerm = x - xCenter, pTerm = p - pCenter;
	ans[0] = sl[j][k][0] + sl_x[j][k][0] * xTerm + sl_p[j][k][0] * pTerm;
	ans[1] = sl[j][k][1] + sl_x[j][k][1] * xTerm + sl_p[j][k][1] * pTerm;
}

/*
 * Helper function for calculating fluxes: function f
 */
double fFcn_Helper(double uFcnVal, double x) {
	return uFcnVal * x;
}

/*
 * Helper function for calculating fluxes: function g
 */
double gFcn_Helper(double omegaFcnVal, double x) {
	return omegaFcnVal * x;
}

/*
 * Calculate the fluxes
 */
void calcFluxesCU() {
	rcstrSolnCU();
	for (int j = 0; j < numCellsX; j++)
		for (int k = 0; k < numCellsP; k++) {
			// The center of this cell
			double xCenter = getCellCenterX(j, k), pCenter = getCellCenterP(j, k);
			// Grid coordinates
			double xRight = getCellRightX(j, k), pTop = getCellTopP(j, k);

			/* Reconstruct solution values at the interfaces, i.e. calculate u_{j, k}^N,
			   u_{j, k + 1}^S, u_{j, k}^E, u_{j + 1, k}^W*/
			double slN[2], slSNext[2], slE[2], slWNext[2];
			calcRcSolnCU(slN, j, k, xCenter, pTop);
			calcRcSolnCU(slSNext, j, k + 1, xCenter, pTop);
			calcRcSolnCU(slE, j, k, xRight, pCenter);
			calcRcSolnCU(slWNext, j + 1, k, xRight, pCenter);

			// u and omega values
			double uFcnValRight = u_fcn(xRight, pCenter),
					omegaFcnValTop = omega_fcn(xCenter, pTop);

			// Calculate the speeds of the interfaces
			double aPlus, aMinus, bPlus, bMinus;
			// Calculate aPlus and aMinus
			double eVal = uFcnValRight;
			aPlus = max(eVal, 0);
			aMinus = min(eVal, 0);
			// Calculate bPlus and bMinus
			eVal = omegaFcnValTop;
			bPlus = max(eVal, 0);
			bMinus = min(eVal, 0);

			// Calculate the fluxes
			if (aPlus != aMinus) {
				for (int ii = 0; ii < 2; ii++)
					Hx[j][k][ii] = (aPlus * fFcn_Helper(uFcnValRight, slE[ii]) -
							aMinus * fFcn_Helper(uFcnValRight, slWNext[ii]) +
							aPlus * aMinus * (slWNext[ii] - slE[ii])) / (aPlus - aMinus);
			} else {
				Hx[j][k][0] = 0;
				Hx[j][k][1] = 0;
			}
			if (bPlus != bMinus) {
				for (int ii = 0; ii < 2; ii++)
					Hp[j][k][ii] = (bPlus * gFcn_Helper(omegaFcnValTop, slN[ii]) -
							bMinus * gFcn_Helper(omegaFcnValTop, slSNext[ii]) +
							bPlus * bMinus * (slSNext[ii] - slN[ii])) / (bPlus - bMinus);
			} else {
				Hp[j][k][0] = 0;
				Hp[j][k][1] = 0;
			}

			/*if (bPlus - bMinus == 0 || aPlus - aMinus == 0) {
				printf("j = %d, k = %d\n", j, k);
				printf("aPlus = %f, aMinus= %f,\n bPlus = %f, bMinus = %f\n", aPlus, aMinus, bPlus, bMinus);
			}*/
		}
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Runge-Kutta Method: RHS Calculation
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void calcRHS_RK_CU(double ans[2], int j, int k) {
	for (int ii = 0; ii < 2; ii++) {
		ans[0] += - (Hx[j][k][ii] + Hx[j - 1][k][ii]) * DxInv
				- (Hp[j][k][ii] - Hp[j][k - 1][ii]) * DpInv;
		ans[1] += 0;
	}
}

#endif /* CUFLUXES_H_ */
