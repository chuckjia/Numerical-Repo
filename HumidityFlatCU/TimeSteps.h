/*
 * TimeSteps.h
 *
 *  Created on: Jun 20, 2017
 *      Author: chuckjia
 */

#ifndef TIMESTEPS_H_
#define TIMESTEPS_H_
#include "CUFluxes.h"
#include "GodunovFluxes.h"

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Select Numerical Scheme
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void (*flux_Scheme)();
void (*RHS_RK_Scheme)(double ans[2], int j, int k);

void selectScheme() {
	if (numericalScheme == 1) {
		flux_Scheme = &calcFluxesCU;
		RHS_RK_Scheme = &calcRHS_RK_CU;
	} else {
		flux_Scheme = &calcFluxesGodunov;
		RHS_RK_Scheme = &addRHS_RK_Godunov;
	}
}

/*
 * Calculate Fluxes: Wrapper Function
 */
void calcFluxes() {
	(*flux_Scheme)();
}

/*
 * Runge-Kutta Method RHS: Wrapper Function
 */
void addRHS_RK(double ans[2], int j, int k) {
	(*RHS_RK_Scheme)(ans, j, k);
}

void forwardEuler() {
	// If using Godunov, we need to cache the u and omega values on the sides
	if (numericalScheme == 0)
		fillCache_uomegaSideVals_Godunov();

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
		for (int j = 1; j < lastIndexX; j++)
			for (int k = 1; k < lastIndexP; k++) {
				// Current solution
				double T = sl[j][k][0], q = sl[j][k][1];
				double x = getCellCenterX(j, k), p = getCellCenterP(j, k);
				// Initiate the RHS of the RK formula
				double RHS_RK[2] = {0, 0};
				// Add fluxes
				addRHS_RK(RHS_RK, j, k);
				// Add the source values to RHS
				addSourceFcn(RHS_RK, T, q, x, p, t, j, k);
				// RK iteration
				for (int ii = 0; ii < 2; ii++)
					sl[j][k][ii] += Dt * RHS_RK[ii];;
			}
		(*enforceBoundaryCondPtr)();
	}
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Runge-Kutta Method: 2nd Order
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void rk2() {
	// If using Godunov, we need to cache the u and omega values on the sides
	if (numericalScheme == 0)
		fillCache_uomegaSideVals_Godunov();

	// Initialize the intermediate solution
	double slCurr[numCellsX][numCellsP][2];

	// Some temporary constants
	double halfDt = 0.5 * Dt;

	// Message on progress
	printf("\n- Using the 2nd order Runge-Kutta method.\n");
	// Temporary value for the progress message
	double timer_factor1_CONST = 100.0 / finalTime;

	for (int tt = 0; tt < numTimeSteps; tt++) {
		double t = Dt * tt;
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
				double RHS_RK[2] = {0, 0};
				// Add fluxes
				addRHS_RK_Godunov(RHS_RK, j, k);
				// Add the source values to RHS
				addSourceFcn(RHS_RK, T, q, x, p, t, j, k);
				//printf("Before = %f\n", sl[j][k][0]);
				// RK iteration
				for (int ii = 0; ii < 2; ii++)
					sl[j][k][ii] += halfDt * RHS_RK[ii];
				//printf("After = %f\n\n", sl[j][k][0]);
			}
		(*enforceBoundaryCondPtr)();

		// Runge-Kutta step 2
		calcFluxes();
		t = t + halfDt;
		for (int j = 1; j < lastIndexX; j++)
			for (int k = 1; k < lastIndexP; k++) {
				// Current solution
				double T = sl[j][k][0], q = sl[j][k][1];
				double x = getCellCenterX(j, k), p = getCellCenterP(j, k);
				// Initiate the RHS of the RK formula
				double RHS_RK[2] = {0, 0};
				// Add fluxes
				addRHS_RK(RHS_RK, j, k);
				// Add the source values to RHS
				addSourceFcn(RHS_RK, T, q, x, p, t, j, k);
				// RK iteration
				for (int ii = 0; ii < 2; ii++)
					sl[j][k][ii] = Dt * RHS_RK[ii] + slCurr[j][k][ii];
			}
		(*enforceBoundaryCondPtr)();
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
		rk2();
	else
		forwardEuler();

	// Measure execution time
	endTime = clock();
	double cpuTimeUsed = (double) (endTime - startTime) * 0.000001 ;
	printf("\r  [Completed]: 100%%. (%1.2fs)\n", cpuTimeUsed);
}

#endif /* TIMESTEPS_H_ */
