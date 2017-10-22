/*
 * TimeSteps.h
 *
 *  Created on: Jun 20, 2017
 *      Author: chuckjia
 */

#ifndef TIMESTEPS_H_
#define TIMESTEPS_H_
#include "SetUpTests.h"

void forwardEuler() {
	// Message on progress
	printf("\n- Using the forward Euler method.\n");
	// Temporary value for the progress message
	double timer_factor1_CONST = 100.0 / finalTime;

	for (int tt = 0; tt < numTimeSteps; tt++) {
		double t = Dt * tt;
		// Printing progress
		printf("\r  [Progress]:%5.1f%%", t * timer_factor1_CONST);
		fflush(stdout);

		// Runge-Kutta step 1
		calcFluxes();
		for (int j = 1; j < lastGhostX; j++)
			for (int k = 1; k < lastGhostP; k++) {
				// Current solution
				double T = soln[j][k][0], q = soln[j][k][1];
				double x = getCellCenterX(j, k), p = getCellCenterP(j, k);
				// Initiate the RHS of the RK formula
				double RHS_RK[2] = {0, 0};
				// Add fluxes
				addRHS_RK(RHS_RK, j, k);
				// Add the source values to RHS
				addSourceFcn(RHS_RK, T, q, x, p, t, j, k);
				// RK iteration
				for (int ii = 0; ii < 2; ii++)
					soln[j][k][ii] += Dt * RHS_RK[ii];;
			}
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
	if (timeScheme == 1)
		forwardEuler();

	// Measure execution time
	endTime = clock();
	double cpuTimeUsed = (double) (endTime - startTime) * 0.000001 ;
	printf("\r  [Completed]: 100%%. (%1.2fs)\n", cpuTimeUsed);
}

#endif /* TIMESTEPS_H_ */
