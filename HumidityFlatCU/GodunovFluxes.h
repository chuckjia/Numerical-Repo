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
/*
 * Cache for u and omega side values: component 0 is u at cell top side;
 * component 1 is omega at cell right side
 */
double uomega_Cache_Godunov[numCellsX][numCellsP][2];

/*
 * Fill in cache for u and omega values at cell sides
 */
void fillCache_uomegaSideVals_Godunov() {
	for (int i = 0; i < numCellsX; i++)
		for (int j = 0; j < numCellsP; j++) {
			 uomega_Cache_Godunov[i][j][0] = Dp * u_fcn(getCellRightX(i, j), getCellCenterP(i, j));
			// uomega_Cache_Godunov[i][j][0] = Dp * u_fcn(getCellCenterX(i, j), getCellCenterP(i, j));
			 uomega_Cache_Godunov[i][j][1] = Dx * omega_fcn(getCellCenterX(i, j), getCellTopP(i, j));
			// uomega_Cache_Godunov[i][j][1] = Dx * omega_fcn(getCellCenterX(i, j), getCellCenterP(i, j));
		}
}

void calcFluxesGodunov() {
	// Calculate G fluxes
	for (int i = 1; i < lastIndexX; i++)
		for (int j = 0; j < lastIndexP; j++) {
			double coeff = uomega_Cache_Godunov[i][j][1];
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
	for (int i = 0; i < lastIndexX; i++)
		for (int j = 1; j < lastIndexP; j++) {
			double coeff = uomega_Cache_Godunov[i][j][0];
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

void addRHS_RK_Godunov(double ans[2], int i, int j) {
	double cellVolVal = getCellVol(i, j);
	for (int ii = 0; ii < 2; ii++) {
		ans[ii] += - (GFlux[i][j][ii] - GFlux[i][j - 1][ii] +
				FFlux[i][j][ii] - FFlux[i - 1][j][ii]) / cellVolVal;
	}
}

#endif /* GODUNOVFLUXES_H_ */
