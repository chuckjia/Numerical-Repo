/*
 * Conditions.h
 *
 *  Created on: Jun 21, 2017
 *      Author: chuckjia
 */

#ifndef CONDITIONS_H_
#define CONDITIONS_H_
#include "Models.h"

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Wrapper Functions for Tests
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

double (*calcExactT)(double x, double p, double t, int j, int k);
double (*sourceFcnPtr)(double T, double q, double x, double p, double t,
		int j, int k);

void setUpTests() {
	if (testNumber == 1) {
		prep_Test1();
		calcExactT = &soln_T_Test1;
		sourceFcnPtr = &source1_test1;
	} else {
		prep_Test2();
		calcExactT = &soln_T_Test2;
		sourceFcnPtr = &source1_test2;
	}
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Initial Conditions
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

/*
 * Initialize the solution
 */
double sl[numCellsX][numCellsP][2];

/*
 * Set initial conditions
 */
void setInitCond() {
	// Measure execution time
	clock_t startTime, endTime;
	startTime = clock();

	// Main body
	for (int j = 0; j < numCellsX; j++)
		for (int k = 0; k < numCellsP; k++) {
			double x = getCellCenterX(j, k), p = getCellCenterP(j, k);
			sl[j][k][0] = (*calcExactT)(x, p, 0, j, k);
			sl[j][k][1] = 0;
		}

	// Measure execution time
	endTime = clock();
	double cpuTimeUsed = (double) (endTime - startTime) * 0.001;
	printf("\n- Initial conditions set. (%1.1fms)\n", cpuTimeUsed);
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Boundary Conditions
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

const int secLastIndexX = lastIndexX - 1;  // Second last index (x-direction)
const int secLastIndexP = lastIndexP - 1;  // Second last index (p-direction)

/*
 * Neumann condition
 */
void neumannCond() {
	// When j = 0 or lastIndexX
	for (int k = 1; k < lastIndexP; k++) {
		// When j = 0
		sl[0][k][0] = sl[1][k][0];
		sl[0][k][1] = sl[1][k][1];
		// When j = lastIndexX
		sl[lastIndexX][k][0] = sl[secLastIndexX][k][0];
		sl[lastIndexX][k][1] = sl[secLastIndexX][k][1];
	}
	// When k = 0 or lastIndexP
	for (int j = 1; j < lastIndexX; j++) {
		// When k = 0
		sl[j][0][0] = sl[j][1][0];
		sl[j][0][1] = sl[j][1][1];
		// When k = lastIndexP
		sl[j][lastIndexP][0] = sl[j][secLastIndexP][0];
		sl[j][lastIndexP][1] = sl[j][secLastIndexP][1];
	}
}

void dirichletCond() {
	// When j = 0 or lastIndexX
	for (int k = 1; k < lastIndexP; k++) {
		// When j = 0
		sl[0][k][0] = 0;
		sl[0][k][1] = 0;
		// When j = lastIndexX
		sl[lastIndexX][k][0] = 0;
		sl[lastIndexX][k][1] = sl[secLastIndexX][k][1];
	}
	// When k = 0 or lastIndexP
	for (int j = 1; j < lastIndexX; j++) {
		// When k = 0
		sl[j][0][0] = 0;
		sl[j][0][1] = 0;
		// When k = lastIndexP
		sl[j][lastIndexP][0] = 0;
		sl[j][lastIndexP][1] = 0;
	}
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Source Terms
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

/*
 * Source function: wrapper
 */
void calcSourceFcn(double ans[2], double T, double q, double x, double p, double t, int j, int k) {
	ans[0] = (*sourceFcnPtr)(T, q, x, p, t, j, k);
	ans[1] = 0;
}



#endif /* CONDITIONS_H_ */
