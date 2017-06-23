/*
 * Fluxes.h
 *
 *  Created on: Jun 22, 2017
 *      Author: chuckjia
 */

#ifndef FLUXES_H_
#define FLUXES_H_
#include "Conditions.h"

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Initialize the Fluxes
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Initialize the fluxes
double Hx[numCellsX][numCellsP][2];
double Hp[numCellsX][numCellsP][2];

// Initialize parameters for the reconstructed solutions
double uX[numCellsX][numCellsP][2];
double uP[numCellsX][numCellsP][2];

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Central Upwind: Calculate the Fluxes
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

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
			uX[j][k][0] = minmod3(theta_CONST * (slRight1 - slThis1) * DxInv,
					0.5 * (slRight1 - slLeft1) * DxInv,
					theta_CONST * (slThis1 - slLeft1) * DxInv);
			uX[j][k][1] = minmod3(theta_CONST * (slRight2 - slThis2) * DxInv,
					0.5 * (slRight2 - slLeft2) * DxInv,
					theta_CONST * (slThis2 - slLeft2) * DxInv);
			uP[j][k][0] = minmod3(theta_CONST * (slUpper1 - slThis1) * DpInv,
					0.5 * (slUpper1 - slLower1) * DpInv,
					theta_CONST * (slThis1 - slLower1) * DpInv);
			uP[j][k][1] = minmod3(theta_CONST * (slUpper2 - slThis2) * DpInv,
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
			uP[j][k][0] = minmod3(theta_CONST * (slUpper1 - slThis1) * DpInv,
					0.5 * (slUpper1 - slLower1) * DpInv,
					theta_CONST * (slThis1 - slLower1) * DpInv);
			uP[j][k][1] = minmod3(theta_CONST * (slUpper2 - slThis2) * DpInv,
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
			uX[j][k][0] = minmod3(theta_CONST * (slRight1 - slThis1) * DxInv,
					0.5 * (slRight1 - slLeft1) * DxInv,
					theta_CONST * (slThis1 - slLeft1) * DxInv);
			uX[j][k][1] = minmod3(theta_CONST * (slRight2 - slThis2) * DxInv,
					0.5 * (slRight2 - slLeft2) * DxInv,
					theta_CONST * (slThis2 - slLeft2) * DxInv);
		}
	}
}

/*
 * Calculate the reconstructed solution value
 */
void calcRcstrSolnCU(double ans[2], int j, int k, double x, double p) {
	double xCenter = getCenterX(j, k), pCenter = getCenterP(j, k);
	double xTerm = x - xCenter, pTerm = p - pCenter;
	ans[0] = sl[j][k][0] + uX[j][k][0] * xTerm + uP[j][k][0] * pTerm;
	ans[1] = sl[j][k][1] + uX[j][k][1] * xTerm + uP[j][k][1] * pTerm;
}

/*
 * Calculate the fluxes
 */
void calcFluxesCU() {
	rcstrSolnCU();
	for (int j = 0; j < numCellsX; j++)
		for (int k = 0; k < numCellsP; k++) {
			// The center of this cell
			double xCenter = getCenterX(j, k), pCenter = getCenterP(j, k);

		}
}

#endif /* FLUXES_H_ */



























