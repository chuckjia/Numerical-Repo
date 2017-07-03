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
 * Initialize the solution
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

double sl[numCellsX][numCellsP][2];

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Wrapper Functions for Tests
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

double (*initTFcnPtr)(double x, double p, double t, int j, int k);
double (*initqFcnPtr)(double x, double p, double t, int j, int k);
void (*sourceFcnPtr)(double ans[2], double T, double q, double x, double p, double t,
		int j, int k);
void (*addBoundaryCondPtr)();

/*
 * Pre-declaration for the boundary conditions
 */
void neumannCond();
void dirichletCond();

void setUpTests() {
	if (testNumber == 0) {
		initTFcnPtr = &initTOrig;
		initqFcnPtr = &initqOrig;
		sourceFcnPtr = &source_Orig;
		addBoundaryCondPtr = &neumannCond;
	} else if (testNumber == 1) {
		prep_Test1();
		initTFcnPtr = &soln_T_Test1;
		initqFcnPtr = &zeroInit;
		sourceFcnPtr = &source_Test1;
		addBoundaryCondPtr = &neumannCond;
	} else if (testNumber == 2) {
		prep_Test2();
		initTFcnPtr = &soln_T_Test2;
		initqFcnPtr = &zeroInit;
		sourceFcnPtr = &source_Test2;
		addBoundaryCondPtr = &dirichletCond;
	} else {
		prep_Test3();
		initTFcnPtr = &soln_T_Test3;
		initqFcnPtr = &zeroInit;
		sourceFcnPtr = &source_Test3;
		addBoundaryCondPtr = &dirichletCond;
	}
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

/*
 * Dirichlet condition
 */
void dirichletCond() {
	// When j = 0 or lastIndexX
	for (int k = 1; k < lastIndexP; k++) {
		// When j = 0
		sl[0][k][0] = 0;
		sl[0][k][1] = 0;
		// When j = lastIndexX
		sl[lastIndexX][k][0] = 0;
		sl[lastIndexX][k][1] = 0;
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
 * Initial Conditions
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

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
			sl[j][k][0] = (*initTFcnPtr)(x, p, 0, j, k);
			sl[j][k][1] = (*initqFcnPtr)(x, p, 0, j, k);
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
void calcSourceFcn(double ans[2], double T, double q, double x, double p, double t, int j, int k) {
	(*sourceFcnPtr)(ans, T, q, x, p, t, j, k);
}



#endif /* CONDITIONS_H_ */
