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

void (*flux_Scheme)(double sl[numCellsX][numCellsP][2]);
void (*RHS_RK_Scheme)(double ans[2], int j, int k);

void selectScheme() {
	if (numericalScheme == 1) {
		flux_Scheme = &calcFluxesCU;
		RHS_RK_Scheme = &calcRHS_RK_CU;
	} else {
		flux_Scheme = &calcFluxesGodunov;
		RHS_RK_Scheme = &calcRHS_RK_Godunov;
	}
}

/*
 * Calculate Fluxes: Wrapper Function
 */
void calcFluxes(double sl[numCellsX][numCellsP][2]) {
	(*flux_Scheme)(sl);
}

/*
 * Runge-Kutta Method RHS: Wrapper Function
 */
void calcRHS_RK(double ans[2], int j, int k) {
	(*RHS_RK_Scheme)(ans, j, k);
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Runge-Kutta Method: 2nd Order
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

const double timer_factor1_CONST = 100 / finalTime;

void rk2() {
	// Initialize the intermediate solution
	double slCurr[numCellsX][numCellsP][2];
	// Some temporary constansts
	double halfDt = 0.5 * Dt;
	// Message on progress
	printf("\n- Using the 2nd order Runge-Kutta method.\n");
	for (int ii = 1; ii <= numTimeSteps; ii++) {
		double t = Dt * ii;
		// Printing progress
		printf("\r  [Progress]:%5.1f%%", t * timer_factor1_CONST);
		fflush(stdout);

		// Runge-Kutta step 1
		calcFluxes(sl);
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
					sl[j][k][ii] += halfDt * RHS_RK[ii];
			}
		neumannCond();

		// Runge-Kutta step 2
		calcFluxes(sl);
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
 * Runge-Kutta Method: 2nd Order Fractional Steps
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void rk2Fractional() {
	// Initialize the intermediate solution
	double slCurr[numCellsX][numCellsP][2];
	// Some temporary constansts
	double halfDt = 0.5 * Dt;
	// Message on progress
	printf("\n- Using the 2nd order Runge-Kutta method.\n");
	for (int tt = 1; tt <= numTimeSteps; tt++) {
		double t = Dt * tt;
		// Printing progress
		printf("\r  [Progress]:%5.1f%%", t * timer_factor1_CONST);
		fflush(stdout);

		if (tt % 2 == 1) {
			// Runge-Kutta step 1
			calcFluxes(sl);
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
					calcRHS_RK(RHS_RK, j, k);
					// RK iteration
					for (int ii = 0; ii < 2; ii++)
						sl[j][k][ii] += halfDt * RHS_RK[ii];
				}
			neumannCond();

			// Runge-Kutta step 2
			calcFluxes(sl);
			for (int j = 1; j < numCellsX; j++)
				for (int k = 1; k < numCellsP; k++) {
					// Current solution
					double T = sl[j][k][0], q = sl[j][k][1];
					double x = getCellCenterX(j, k), p = getCellCenterP(j, k);
					// Initiate the RHS of the RK formula
					double RHS_RK[2] = {0, 0};
					// Add fluxes
					calcRHS_RK(RHS_RK, j, k);
					// RK iteration
					for (int ii = 0; ii < 2; ii++)
						sl[j][k][ii] = slCurr[j][k][ii] + Dt * RHS_RK[ii];
				}
			neumannCond();
		} else {
			for (int j = 1; j < lastIndexX; j++)
				for (int k = 1; k < lastIndexP; k++) {
					// Current solution
					double T = sl[j][k][0], q = sl[j][k][1];
					double x = getCellCenterX(j, k), p = getCellCenterP(j, k);
					// Initiate the RHS of the RK formula
					double RHS_RK[2] = {0, 0};
					// Add fluxes
					calcSourceFcn(RHS_RK, T, q, x, p, t, j, k);
					// RK iteration
					for (int ii = 0; ii < 2; ii++)
						sl[j][k][ii] += Dt * RHS_RK[ii];
				}
			neumannCond();
		}
	}
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Runge-Kutta Method: 4th Order
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void rk4() {
	// Initialize the intermediate solution
	double slCurr[numCellsX][numCellsP][2];
	double kVec[numCellsX][numCellsP][2];
	// Some temporary constants
	double halfDt = 0.5 * Dt;
	double DtOver6 = Dt / 6;

	// Message on progress
	printf("\n- Using the 4th order Runge-Kutta method.\n");
	for (int tt = 1; tt <= numTimeSteps; tt++) {
		double t = Dt * tt;
		// Printing progress
		printf("\r  [Progress]:%5.1f%%", t * timer_factor1_CONST);
		fflush(stdout);

		// Runge-Kutta step 1
		calcFluxes(sl);
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
				for (int ii = 0; ii < 2; ii++) {
					kVec[j][k][ii] = RHS_RK[ii];  // Adding k1
					sl[j][k][ii] += halfDt * RHS_RK[ii];  // To be used in Step 2
				}
			}
		neumannCond();

		// Runge-Kutta step 2
		double tCurr = t + halfDt;
		calcFluxes(sl);
		for (int j = 1; j < lastIndexX; j++)
			for (int k = 1; k < lastIndexP; k++) {
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
				for (int ii = 0; ii < 2; ii++) {
					kVec[j][k][ii] += 2 * RHS_RK[ii];  // Adding k2
					sl[j][k][ii] = slCurr[j][k][ii] + halfDt * RHS_RK[ii];  // To be used in Step 3
				}
			}
		neumannCond();

		// Runge-Kutta step 3
		calcFluxes(sl);
		for (int j = 1; j < lastIndexX; j++)
			for (int k = 1; k < lastIndexP; k++) {
				// Current solution
				double T = sl[j][k][0], q = sl[j][k][1];
				double x = getCellCenterX(j, k), p = getCellCenterP(j, k);
				// Initiate the RHS of the RK formula
				double RHS_RK[2];
				// Add the source values to RHS
				calcSourceFcn(RHS_RK, T, q, x, p, tCurr, j, k);
				// Add fluxes
				calcRHS_RK(RHS_RK, j, k);
				// RK iteration
				for (int ii = 0; ii < 2; ii++) {
					kVec[j][k][ii] += 2 * RHS_RK[ii];  // Adding k3
					sl[j][k][ii] = slCurr[j][k][ii] + RHS_RK[ii];  // To be used in Step 4
				}
			}
		neumannCond();

		// Runge-Kutta step 4
		tCurr = t + Dt;
		calcFluxes(sl);
		for (int j = 1; j < lastIndexX; j++)
			for (int k = 1; k < lastIndexP; k++) {
				// Current solution
				double T = sl[j][k][0], q = sl[j][k][1];
				double x = getCellCenterX(j, k), p = getCellCenterP(j, k);
				// Initiate the RHS of the RK formula
				double RHS_RK[2];
				// Add the source values to RHS
				calcSourceFcn(RHS_RK, T, q, x, p, tCurr, j, k);
				// Add fluxes
				calcRHS_RK(RHS_RK, j, k);
				// RK iteration
				for (int ii = 0; ii < 2; ii++) {
					kVec[j][k][ii] += RHS_RK[ii];  // Adding k4
					sl[j][k][ii] = slCurr[j][k][ii] + DtOver6 * kVec[j][k][ii];
				}
			}
		neumannCond();
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
		rk4();

	// Measure execution time
	endTime = clock();
	double cpuTimeUsed = (double) (endTime - startTime) * 0.000001 ;
	printf("\r  [Completed]: 100%%. (%1.2fs)\n", cpuTimeUsed);
}

#endif /* TIMESTEPS_H_ */
