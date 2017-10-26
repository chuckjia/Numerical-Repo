/*
 * Models.h
 *
 *  Created on: Oct 8, 2017
 *      Author: chuckjia
 */

#ifndef MODELS_H_
#define MODELS_H_
#include "Constants.h"

/* ===== ===== ===== ===== ===== =====
 * Model 1
 * ===== ===== ===== ===== ===== ===== */

const double c1_Test1 = M_PI / 5;

// Initial Condition
double exactSoln_Test1(double x, double t) {
	return sin(c1_Test1 * (x - a * t));
}

/* ===== ===== ===== ===== ===== =====
 * Model 2
 * ===== ===== ===== ===== ===== ===== */

const double a_Test2 = 1;  // Velocity coefficient a
double domainCenter = (xf + x0) / 2;

// Initial Condition
double exactSoln_Test2(double x, double t) {
	if (x < domainCenter + a_Test2 * t)
		return 1;
	return 0;
}

void enforceBC_Test2() {
	u_sl[0] = 1;
	u_sl[lastGhostIndex] = u_sl[lastRealIndex];
}

/* ===== ===== ===== ===== ===== =====
 * Initial Conditions
 * ===== ===== ===== ===== ===== ===== */

double (*exactSolnPtr)(double x, double t);

void enforceIC() {
	for (int i = 0; i < numCells; i++) {
		double x = getCellCenter(i);
		u_sl[i] = (*exactSolnPtr)(x, 0);
	}
}

/* ===== ===== ===== ===== ===== =====
 * Boundary Conditions
 * ===== ===== ===== ===== ===== ===== */

void (*enforceBCPtr)();

void enforcePeriodicBC() {
	u_sl[0] = u_sl[lastRealIndex];
	u_sl[lastGhostIndex] = u_sl[1];
}

/* ===== ===== ===== ===== ===== =====
 * Assign Parameters
 * ===== ===== ===== ===== ===== ===== */

void selectModel() {
	// if (modelNum == 1)
	exactSolnPtr = &exactSoln_Test1;
	enforceBCPtr = &enforcePeriodicBC;
	if (modelNum == 2) {
		exactSolnPtr = &exactSoln_Test2;
		enforceBCPtr = &enforceBC_Test2;
	}
}


#endif /* MODELS_H_ */
