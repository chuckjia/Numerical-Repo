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
 * Wrapper Functions for Model Selections
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

double (*initTFcnPtr)(double x, double p, double t, int j, int k);
double (*initqFcnPtr)(double x, double p, double t, int j, int k);
void (*sourceFcnPtr)(double ans[2], double T, double q, double x, double p, double t,
		int j, int k);
void (*enforceBoundaryCondPtr)();

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Flux Functions
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

/*
 * The flux function f
 */
double f_fcn(double input, double x, double p) {
	return u_fcn(x, p) * input;
}

/*
 * The flux function g
 */
double g_fcn(double input, double x, double p) {
	return omega_fcn(x, p) * input;
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
	for (int k = 0; k < numCellsP; k++) {
		// When j = 0
		sl[0][k][0] = sl[1][k][0];
		sl[0][k][1] = sl[1][k][1];
		// When j = lastIndexX
		sl[lastIndexX][k][0] = sl[secLastIndexX][k][0];
		sl[lastIndexX][k][1] = sl[secLastIndexX][k][1];
	}
	// When k = 0 or lastIndexP
	for (int j = 0; j < numCellsX; j++) {
		// When k = 0
		sl[j][0][0] = sl[j][1][0];
		sl[j][0][1] = sl[j][1][1];
		// When k = lastIndexP
		sl[j][lastIndexP][0] = sl[j][secLastIndexP][0];
		sl[j][lastIndexP][1] = sl[j][secLastIndexP][1];
	}
}

/*
 * Neumann condition
 */
void neumannCondLeftRight() {
	// When j = 0 or lastIndexX
	for (int k = 0; k < numCellsP; k++) {
		// When j = 0
		sl[0][k][0] = sl[1][k][0];
		sl[0][k][1] = sl[1][k][1];
		// When j = lastIndexX
		sl[lastIndexX][k][0] = sl[secLastIndexX][k][0];
		sl[lastIndexX][k][1] = sl[secLastIndexX][k][1];
	}
}

/*
 * Dirichlet condition
 */

// Dirichlet condition constant
double boundaryVal = 0;

// Enforce the Dirichlet condition
void dirichletCond() {
	// When j = 0 or lastIndexX
	for (int k = 0; k < numCellsP; k++) {
		// When j = 0
		sl[0][k][0] = boundaryVal;
		sl[0][k][1] = boundaryVal;
		// When j = lastIndexX
		sl[lastIndexX][k][0] = boundaryVal;
		sl[lastIndexX][k][1] = boundaryVal;
	}
	// When k = 0 or lastIndexP
	for (int j = 0; j < numCellsX; j++) {
		// When k = 0
		sl[j][0][0] = boundaryVal;
		sl[j][0][1] = boundaryVal;
		// When k = lastIndexP
		sl[j][lastIndexP][0] = boundaryVal;
		sl[j][lastIndexP][1] = boundaryVal;
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
void addSourceFcn(double ans[2], double T, double q, double x, double p, double t, int j, int k) {
	(*sourceFcnPtr)(ans, T, q, x, p, t, j, k);
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Set Up The Test Parameters
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void setUpTests() {
	if (testNumber == 0) {
		prep_Orig();
		initTFcnPtr = &initTOrig;
		initqFcnPtr = &initqOrig;
		sourceFcnPtr = &source_Orig;
		enforceBoundaryCondPtr = &neumannCond;
	} else if (testNumber == 1) {
		prep_Test1();
		initTFcnPtr = &soln_T_Test1;
		initqFcnPtr = &zeroInit;
		sourceFcnPtr = &source_Test1;
		enforceBoundaryCondPtr = &dirichletCond;
	} else if (testNumber == 2) {
		prep_Test2();
		initTFcnPtr = &soln_T_Test2;
		initqFcnPtr = &zeroInit;
		sourceFcnPtr = &source_Test2;
		boundaryVal = 0;
		enforceBoundaryCondPtr = &dirichletCond;
	} else if (testNumber == 3) {
		prep_Test3();
		initTFcnPtr = &soln_T_Test3;
		initqFcnPtr = &soln_q_Test3;
		sourceFcnPtr = &source_Test3;
		boundaryVal = 0;
		enforceBoundaryCondPtr = &dirichletCond;
	} else if (testNumber == 4) {
		prep_Test4();
		initTFcnPtr = &soln_T_Test4;
		initqFcnPtr = &soln_q_Test4;
		sourceFcnPtr = &source_Test4;
		boundaryVal = 0;
		enforceBoundaryCondPtr = &neumannCond;
	} else if (testNumber == 5) {
		prep_Test5();
		initTFcnPtr = &soln_T_Test5;
		initqFcnPtr = &soln_q_Test5;
		sourceFcnPtr = &source_Test5;
		boundaryVal = 0;
		enforceBoundaryCondPtr = &neumannCondLeftRight;
	}
}



#endif /* CONDITIONS_H_ */
