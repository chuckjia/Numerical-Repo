/*
 * TimeSteps.h
 *
 *  Created on: Jun 20, 2017
 *      Author: chuckjia
 */

#ifndef TIMESTEPS_H_
#define TIMESTEPS_H_
#include "CUFluxes.h"
#include "GodnuvFluxes.h"

/*
 * Initialize the intermediate solution
 */
double slCurr[numCellsX][numCellsP][2];

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Select Numerical Scheme
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

/*
 * Calculate Fluxes: Wrapper Function
 */
void calcFluxes() {
	calcFluxesGodnuv();
}

/*
 * Runge-Kutta Method RHS: Wrapper Function
 */
void calcRHS_RK(double ans[2], int j, int k) {
	calcRHS_RK_Godnuv(ans, j, k);
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Runge-Kutta Method: 2nd Order
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

const double timer_factor1_CONST = 100 / finalTime;

void rk2() {
	// Message on progress
	printf("\n- Using the 2nd order Runge-Kutta method.\n");
	for (int ii = 1; ii <= numTimeSteps; ii++) {
		double t = Dt * ii;
		// Printing progress
		printf("\r  [Progress]:%5.1f%%", t * timer_factor1_CONST);
		fflush(stdout);

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
				double RHS_RK[2];
				// Add the source values to RHS
				calcSourceFcn(RHS_RK, T, q, x, p, t, j, k);
				// Add fluxes
				calcRHS_RK(RHS_RK, j, k);
				// RK iteration
				for (int ii = 0; ii < 2; ii++)
					sl[j][k][ii] += 0.5 * Dt * RHS_RK[ii];
			}
		neumannCond();

		// Runge-Kutta step 2
		calcFluxes();
		for (int j = 1; j < numCellsX; j++)
			for (int k = 1; k < numCellsP; k++) {
				// Current solution
				double T = sl[j][k][0], q = sl[j][k][1];
				double x = getCellCenterX(j, k), p = getCellCenterP(j, k);
				// Initiate the RHS of the RK formula
				double RHS_RK[2];
				// Add the source values to RHS
				calcSourceFcn(RHS_RK, T, q, x, p, t, j, k);
				// Add fluxes
				calcRHS_RK(RHS_RK, j, k);
				// RK iteration
				for (int ii = 0; ii < 2; ii++)
					sl[j][k][ii] = slCurr[j][k][ii] + Dt * RHS_RK[ii];
			}
		neumannCond();
	}
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Runge-Kutta Method: 4th Order
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void rk4() {

}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Time Method: Wrapper
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void timeMethod() {
	// Measure execution time
	clock_t startTime, endTime;
	startTime = clock();

	// Main body
	rk2();

	// Measure execution time
	endTime = clock();
	double cpuTimeUsed = (double) (endTime - startTime) * 0.000001 ;
	printf("\r  [Completed]: 100%%. (%1.2fs)\n", cpuTimeUsed);
}

#endif /* TIMESTEPS_H_ */
