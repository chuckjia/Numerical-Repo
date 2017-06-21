/*
 * TimeMethod.h
 *
 *  Created on: Jun 10, 2017
 *      Author: chuckjia
 */

#ifndef TIMEMETHOD_H_
#define TIMEMETHOD_H_
#include "Fluxes.h"

// Temporarily store the solution for use in the Runge-Kutta method
double solnCurr[numCellsXDir][numCellsPDir][2];

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Calculate the vector R
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void calcRVal(double RVal[2], int j, int k, double t, double DpVal) {
	double sourceVal[2];
	sourceFcn(sourceVal, soln[j][k][0], soln[j][k][1],
			getCenterX(j, k), getCenterP(j, k), t);
	RVal[0] = - (Hx[j][k][0] - Hx[j - 1][k][0]) / Dx -
			(Hp[j][k][0] - Hp[j][k - 1][0]) / DpVal + sourceVal[0];
	RVal[1] = - (Hx[j][k][1] - Hx[j - 1][k][1]) / Dx -
			(Hp[j][k][1] - Hp[j][k - 1][1]) / DpVal + sourceVal[1];
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * 2nd order Runge-Kutta method
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void rk2() {
	for (int ii = 1; ii <= numTimeSteps; ii++){
		double t = ii * Dt;
		printf("%5.1f percent\n", t / finalTime * 100);
		// Step 1
		calcFluxes();
		for (int j = 1; j < numGridPtsXDir; j++) {
			double DpVal = getDp(j);
			for (int k = 1; k < numGridPtsPDir; k++) {
				// Preserve the current solution
				solnCurr[j][k][0] = soln[j][k][0];
				solnCurr[j][k][1] = soln[j][k][1];
				// Calculate R
				double RVal[2];
				calcRVal(RVal, j, k, t, DpVal);
				// Update solution to the intermediate solution
				soln[j][k][0] += 0.5 * Dt * RVal[0];
				soln[j][k][1] += 0.5 * Dt * RVal[1];
			}
		}

		// Step 2
		calcFluxes();
		for (int j = 1; j < numGridPtsXDir; j++) {
			double DpVal = getDp(j);
			for (int k = 1; k < numGridPtsPDir; k++) {
				// Calculate R
				double RVal[2];
				calcRVal(RVal, j, k, t, DpVal);
				// Update the solution for the next time step
				soln[j][k][0] = solnCurr[j][k][0] + Dt * RVal[0];
				soln[j][k][1] = solnCurr[j][k][1] + Dt * RVal[1];
			}
		}

		// The boundary conditions
		// When j = 0 (left side)
		for (int k = 1; k < numGridPtsPDir; k++) {
			soln[0][k][0] = soln[1][k][0];
			soln[0][k][1] = soln[1][k][1];
		}
		// When j = numCellsXDir - 1 (right side)
		int jVal1 = numCellsXDir - 1, jVal2 = numCellsXDir - 2;
		for (int k = 1; k < numGridPtsPDir; k++) {
			soln[jVal1][k][0] = soln[jVal2][k][0];
			soln[jVal1][k][1] = soln[jVal2][k][1];
		}
		// When k = 0 (bottom side)
		for (int j = 1; j < numGridPtsXDir; j++) {
			soln[j][0][0] = soln[j][1][0];
			soln[j][0][1] = soln[j][1][1];
		}
		// When k = numCellsPDir - 1 (top side)
		int kVal1 = numCellsPDir - 1, kVal2 = numCellsPDir - 2;
		for (int j = 1; j < numGridPtsPDir; j++) {
			soln[j][kVal1][0] = soln[j][kVal2][0];
			soln[j][kVal1][1] = soln[j][kVal2][1];
		}
	}

}

#endif /* TIMEMETHOD_H_ */
