/*
 * Models.h
 *
 *  Created on: Sep 3, 2017
 *      Author: chuckjia
 */

#ifndef MODELS_H_
#define MODELS_H_

#include "Constants.h"

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Model Functions
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

/* ----- ----- -----
 * Model 1
 * ----- ----- ----- */

// Coefficients c0, c1, c2 and c3 used in the manufactured solution
const double numOscil_Test1 = 2;  // Number of oscillations in the exact solution
const double c0_Test1 = numOscil_Test1 / (xf - x0);
const double c1_Test1 = c0_Test1 * 2 * M_PI;
const double c2_Test1 = c0_Test1 * 2 * M_PI;
const double c3_Test1 = 2 * M_PI;
const double multiFactor_Test1 = 1;

// The manufactured solutions
double exact_T_Test1(double x, double p, double t) {
	return multiFactor_Test1 * sin(c1_Test1 * x) * sin(c2_Test1 * p) * cos(c3_Test1 * t);
}

double exact_q_Test1(double x, double p, double t) {
	return 0;
}

// Source terms
void addSource_Test1(double ans[2], double T, double q, double x, double p, double t) {
	double xInput = c1_Test1 * x, pInput = c2_Test1 * p, tInput = c3_Test1 * t;
	double xPart = sin(xInput), pPart = sin(pInput), tPart = cos(tInput);
	ans[0] += T * (u_xDer_fcn(x, p) + w_pDer_fcn(x, p)) +
			multiFactor_Test1 * (c1_Test1 * cos(c1_Test1 * x) * pPart * tPart * u_fcn(x, p) +
					xPart * c2_Test1 * cos(c2_Test1 * p) * tPart * w_fcn(x, p) -
					xPart * pPart * c3_Test1 * sin(tInput));
}

/* ----- ----- -----
 * Model 2
 * ----- ----- ----- */

// Coefficients c0, c1, c2 and c3 used in the manufactured solution
const double numOscil_Test2 = 2;
const double c0_Test2 = numOscil_Test2 / (xf - x0);
const double c1_Test2 = c0_Test2 * 2 * M_PI;
const double c2_Test2 = c0_Test2 * 2 * M_PI;
const double c3_Test2 = 2 * M_PI;

// The manufactured solution
double exact_T_Test2(double x, double p, double t) {
	return cos(c1_Test2 * x) * cos(c2_Test2 * p) * cos(c3_Test2 * t);
}

double exact_q_Test2(double x, double p, double t) {
	return 0;
}

// Source terms
void addSource_Test2(double ans[2], double T, double q, double x, double p, double t) {
	double xInput = c1_Test2 * x, pInput = c2_Test2 * p, tInput = c3_Test2 * t;
	double xPart = cos(xInput), pPart = cos(pInput), tPart = cos(tInput);
	ans[0] += -xPart * pPart * c3_Test2 * sin(tInput) +
			T * (u_xDer_fcn(x, p) + w_pDer_fcn(x, p)) -
			c1_Test2 * sin(c1_Test2 * x) * pPart * tPart * u_fcn(x, p) -
			xPart * c2_Test2 * sin(c2_Test2 * p) * tPart * w_fcn(x, p);
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Initial Conditions
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Allocate the space for solutions
double sl[numCellsX][numCellsP][2];
// The initial condition
double (*initTPtr)(double x, double p, double t),
		(*initqPtr)(double x, double p, double t);

// Enforce the initial conditions
void enforceInitCond() {
	// Select the initial conditions
	for (int i = 0; i < numCellsX; i++)
		for (int j = 0; j < numCellsP; j++) {
			double x = getCellCenterX(i, j), p = getCellCenterP(i, j);
			sl[i][j][0] = (*initTPtr)(x, p, 0);
			sl[i][j][1] = (*initqPtr)(x, p, 0);
		}
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Boundary Conditions
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void (*enforceBCPtr)();

// Boundary values for the Dirichlet condition
double boundaryVal = 0;

// Enforce the Dirichlet boundary conditions by extrapolation on Ghost cells
void enforceDirichletBC() {
	for (int kk = 0; kk < 2; kk++) {
		double twiceBoundVal = 2 * boundaryVal;
		// On the left and right boundaries
		for (int j = 1; j < lastGhostIndexP; j++) {
			sl[0][j][kk] = twiceBoundVal - sl[1][j][kk];
			sl[lastGhostIndexX][j][kk] = twiceBoundVal - sl[lastRealIndexX][j][kk];
		}
		// On the top and bottom boundaries
		for (int i = 1; i < lastGhostIndexX; i++) {
			sl[i][0][kk] = twiceBoundVal - sl[i][1][kk];
			sl[i][lastGhostIndexP][kk] = twiceBoundVal - sl[i][lastRealIndexP][kk];
		}
	}
}

// Enforce the Neumann boundary conditions by extrapolation on Ghost cells
void enforceNeumannBC() {
	// On solution T
	for (int j = 1; j < lastGhostIndexP; j++) {
		sl[0][j][0] = sl[1][j][0];
		sl[lastGhostIndexX][j][0] = sl[lastRealIndexX][j][0];
	}
	for (int i = 1; i < lastGhostIndexX; i++) {
		sl[i][0][0] = sl[i][1][0];
		sl[i][lastGhostIndexP][0] = sl[i][lastRealIndexP][0];
	}

	// On solution q
	for (int j = 1; j < lastGhostIndexP; j++) {
		sl[0][j][1] = sl[1][j][1];
		sl[lastGhostIndexX][j][1] = sl[lastRealIndexX][j][1];
	}
	for (int i = 1; i < lastGhostIndexX; i++) {
		sl[i][0][1] = sl[i][1][1];
		sl[i][lastGhostIndexP][1] = sl[i][lastRealIndexP][1];
	}
}

// Enforce the Periodic boundary conditions by extrapolation on Ghost cells
void enforcePeriodicBC() {
	// On solution T
	for (int j = 1; j < lastGhostIndexP; j++) {
		sl[0][j][0] = sl[lastRealIndexX][j][0];
		sl[lastGhostIndexX][j][0] = sl[1][j][0];
	}
	for (int i = 1; i < lastGhostIndexX; i++) {
		sl[i][0][0] = sl[i][lastRealIndexP][0];
		sl[i][lastGhostIndexP][0] = sl[i][1][0];
	}

	// On solution q
	for (int j = 1; j < lastGhostIndexP; j++) {
		sl[0][j][1] = sl[lastRealIndexX][j][1];
		sl[lastGhostIndexX][j][1] = sl[1][j][1];
	}
	for (int i = 1; i < lastGhostIndexX; i++) {
		sl[i][0][1] = sl[i][lastRealIndexP][1];
		sl[i][lastGhostIndexP][1] = sl[i][1][1];
	}
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Source Functions
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void (*addSourceFcnPtr)(double ans[2], double T, double q, double x, double p, double t);

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Model Selection: Assign Values Function Pointers
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void selectModel() {
	if (modelNum == 1) {
		// Initial conditions
		initTPtr = &exact_T_Test1; initqPtr = &exact_q_Test1;
		// Boundary conditions
		boundaryVal = 0; enforceBCPtr = &enforceDirichletBC;
		// Source functions
		addSourceFcnPtr = &addSource_Test1;
	} else if (modelNum == 2) {
		// Initial conditions
		initTPtr = &exact_T_Test2; initqPtr = &exact_q_Test2;
		// Boundary conditions
		enforceBCPtr = &enforceNeumannBC;
		// Source functions
		addSourceFcnPtr = &addSource_Test2;
	}
}

#endif /* MODELS_H_ */
