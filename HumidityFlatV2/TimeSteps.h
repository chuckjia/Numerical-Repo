/*
 * TimeSteps.h
 *
 *  Created on: Jun 20, 2017
 *      Author: chuckjia
 */

#ifndef TIMESTEPS_H_
#define TIMESTEPS_H_
#include "GodunovFluxes.h"

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Wrapper Functions
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void calcFluxes() {
	calcFluxesGodunov();
}

void calcRHS_RK(double ans[2], int i, int j) {
	calcRHS_RK_Godunov(ans, i, j);
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Runge-Kutta Method: 2nd Order
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

const double timer_factor1_CONST = 100 / finalTime;

void rk2_Godunov() {
	// Initialize the values for u and omega values to be used in Godunov
	fillCache_uomegaSideVals_Godunov();
	// Initialize the intermediate solution
	double slCurr[numCellsX][numCellsP][2];
	// Some temporary constants
	double halfDt = 0.5 * Dt;
	// Message on progress
	printf("\n- Using the 2nd order Runge-Kutta method.\n");
	for (int ii = 0; ii < numTimeSteps; ii++) {
		double t = Dt * ii;
		// Printing progress
		/*printf("\r  [Progress]:%5.1f%%", t * timer_factor1_CONST);
		fflush(stdout);*/

		// Runge-Kutta step 1
		calcFluxes();
		for (int j = 1; j < lastIndexX; j++)
			for (int k = 1; k < lastIndexP; k++) {
				// Current solution
				double T = sl[j][k][0], q = sl[j][k][1];
				double x = getCellCenterX(j, k), p = getCellCenterP(j, k);
				// Store the value of the current solution
				slCurr[j][k][0] = T;
				slCurr[j][k][1] = q;
				// Initiate the RHS of the RK formula
				double RHS_RK[2] = {0, 0};
				// Add the source values to RHS
				calcSourceFcn(RHS_RK, T, q, x, p, t, j, k);
				// Add fluxes
				calcRHS_RK(RHS_RK, j, k);
				// RK iteration
				for (int ii = 0; ii < 2; ii++)
					sl[j][k][ii] += halfDt * RHS_RK[ii];
			}
		(*addBoundaryCondPtr)();

		// Runge-Kutta step 2
		calcFluxes();
		for (int j = 1; j < numCellsX; j++)
			for (int k = 1; k < numCellsP; k++) {
				// Current solution
				double T = sl[j][k][0], q = sl[j][k][1];
				double x = getCellCenterX(j, k), p = getCellCenterP(j, k);
				// Initiate the RHS of the RK formula
				double RHS_RK[2] = {0, 0};
				// Add the source values to RHS
				calcSourceFcn(RHS_RK, T, q, x, p, t + halfDt, j, k);
				// Add fluxes
				calcRHS_RK(RHS_RK, j, k);
				// RK iteration
				for (int ii = 0; ii < 2; ii++)
					sl[j][k][ii] = slCurr[j][k][ii] + Dt * RHS_RK[ii];
			}
		(*addBoundaryCondPtr)();
	}
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Time Method: Wrapper
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void timeMethod() {
	// Measure execution time
	clock_t startTime, endTime;
	startTime = clock();

	// Main body
	if (timeScheme == 2)
		rk2_Godunov();
	else
		rk2_Godunov();

	// Measure execution time
	endTime = clock();
	double cpuTimeUsed = (double) (endTime - startTime) * 0.000001 ;
	printf("\r  [Completed]: 100%%. (%1.2fs)\n", cpuTimeUsed);
}

#endif /* TIMESTEPS_H_ */
