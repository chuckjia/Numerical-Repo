/*
 * Conditions.h
 *
 *  Created on: July 10, 2017
 *      Author: chuckjia
 */

#ifndef CONDITIONS_H_
#define CONDITIONS_H_
#include "Models.h"

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Initialize the Solution
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

double soln[numCellsX][numCellsP][2];

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Initial Conditions
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

double (*initTFcnPtr)(double x, double p, double t, int i, int j);
double (*initqFcnPtr)(double x, double p, double t, int i, int j);
void (*addSourceFcnPtr)(double ans[2], double T, double q, double x, double p, double t,
		int i, int j);

/*
 * Set initial conditions
 */
void setInitCond() {
	// Measure execution time
	clock_t startTime, endTime;
	startTime = clock();

	// Main body
	for (int i = 0; i < numCellsX; i++)
		for (int j = 0; j < numCellsP; j++) {
			double x = getCellCenterX(i, j), p = getCellCenterP(i, j);
			soln[i][j][0] = (*initTFcnPtr)(x, p, 0, i, j);
			soln[i][j][1] = (*initqFcnPtr)(x, p, 0, i, j);
		}

	// Measure execution time
	endTime = clock();
	double cpuTimeUsed = (double) (endTime - startTime) * 0.001;
	printf("\n- Initial conditions set. (%1.1fms)\n", cpuTimeUsed);
	// printMatrix2DTimes2_Comp1(numCellsX, numCellsP, sl);
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Source Terms
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

/*
 * Source function: wrapper
 */
void addSourceFcn(double ans[2], double T, double q, double x, double p, double t, int i, int j) {
	(*addSourceFcnPtr)(ans, T, q, x, p, t, i, j);
}

#endif /* CONDITIONS_H_ */
