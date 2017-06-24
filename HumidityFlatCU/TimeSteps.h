/*
 * TimeSteps.h
 *
 *  Created on: Jun 20, 2017
 *      Author: chuckjia
 */

#ifndef TIMESTEPS_H_
#define TIMESTEPS_H_
#include "Fluxes.h"

/*
 * Initialize the intermediate solution
 */
double slCurr[numCellsX][numCellsP][2];

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Runge-Kutta Method: 2nd Order
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void rk2() {
	printf("Using the 2nd order Runge-Kutta method.\n");
	for (int ii = 1; ii <= numTimeSteps; ii++) {
		double t = Dt * ii;
		// Printing progress
		printf("%5.1f per cent\n", t / finalTime * 100);

		// Runge-Kutta step 1
		calcFluxes();
		for (int j = 0; j < numCellsX; j++)
			for (int k = 0; k < numCellsP; k++) {
				// Current solution
				double T = sl[j][k][0], q = sl[j][k][1];
				double xCenter = getCellCenterX(j, k), pCenter = getCellCenterP(j, k);
				// Store the value of the current solution
				slCurr[j][k][0] = T;
				slCurr[j][k][1] = q;
				// Initiate the RHS of the RK formula
				double RHS_RK[2];
				// Add the source values to RHS
				calcSourceFcn(RHS_RK, T, q, xCenter, pCenter, t);
				// Add fluxes
				for (int ii = 0; ii < 2; ii++) {
					RHS_RK[ii] += - (Hx[j][k][ii] + Hx[j - 1][k][ii]) * DxInv
							- (Hp[j][k][ii] - Hp[j][k - 1][ii]) * DpInv;
					sl[j][k][ii] += 0.5 * Dt * RHS_RK[ii];
				}
			}

		// Runge-Kutta step 2
		calcFluxes();
		for (int j = 0; j < numCellsX; j++)
			for (int k = 0; k < numCellsP; k++) {
				// Current solution
				double T = sl[j][k][0], q = sl[j][k][1];
				double xCenter = getCellCenterX(j, k), pCenter = getCellCenterP(j, k);
				// Initiate the RHS of the RK formula
				double RHS_RK[2];
				// Add the source values to RHS
				calcSourceFcn(RHS_RK, T, q, xCenter, pCenter, t);
				// Add fluxes
				for (int ii = 0; ii < 2; ii++) {
					RHS_RK[ii] += - (Hx[j][k][ii] + Hx[j - 1][k][ii]) * DxInv
							- (Hp[j][k][ii] - Hp[j][k - 1][ii]) * DpInv;
					sl[j][k][ii] = slCurr[j][k][ii] + Dt * RHS_RK[ii];
				}
			}
	}
}

#endif /* TIMESTEPS_H_ */
