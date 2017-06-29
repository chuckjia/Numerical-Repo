/*
 * GodunovFluxes.h
 *
 *  Created on: Jun 20, 2017
 *      Author: chuckjia
 */

#ifndef GODUNOVFLUXES_H_
#define GODUNOVFLUXES_H_
#include "Conditions.h"

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Initialize the Fluxes
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Initialize the fluxes
// FFlux[j][k] represents F_{i+1/2,j} and GFlux[j][k] represents G_{i,j+1/2}
double FFlux[numCellsX][numCellsP][2];
double GFlux[numCellsX][numCellsP][2];

void calcFluxesGodunov(double sl[numCellsX][numCellsP][2]) {
	// Calculate G fluxes
	for (int i = 1; i < lastIndexX; i++)
		for (int j = 1; j < lastIndexP; j++) {
			double coeff = Dx * omega_fcn(getCellCenterX(i, j), getCellTopP(i, j));
			double slCheck[2] = {0, 0};
			if (coeff >= 0) {
				slCheck[0] = sl[i][j][0];
				slCheck[1] = sl[i][j][1];
			} else {
				slCheck[0] = sl[i][j + 1][0];
				slCheck[1] = sl[i][j + 1][1];
			}
			GFlux[i][j][0] = coeff * slCheck[0];
			GFlux[i][j][1] = coeff * slCheck[1];
		}
	// Calculate F fluxes
	for (int i = 1; i < lastIndexX; i++)
		for (int j = 1; j < lastIndexP; j++) {
			double coeff = Dp * u_fcn(getCellRightX(i, j), getCellCenterP(i, j));
			double slCheck[2] = {0, 0};
			if (coeff >= 0) {
				slCheck[0] = sl[i][j][0];
				slCheck[1] = sl[i][j][1];
			} else {
				slCheck[0] = sl[i + 1][j][0];
				slCheck[1] = sl[i + 1][j][1];
			}
			FFlux[i][j][0] = coeff * slCheck[0];
			FFlux[i][j][1] = coeff * slCheck[1];
		}
}

void calcRHS_RK_Godunov(double ans[2], int i, int j) {
	for (int ii = 0; ii < 2; ii++) {
		ans[ii] += -1 / cellVol * (GFlux[i][j][ii] - GFlux[i][j - 1][ii] +
				FFlux[i][j][ii] - FFlux[i - 1][j][ii]);
	}
}

#endif /* GODUNOVFLUXES_H_ */
