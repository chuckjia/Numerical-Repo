/*
 * SetUpTests.h
 *
 *  Created on: Jul 14, 2017
 *      Author: chuckjia
 */

#ifndef SETUPTESTS_H_
#define SETUPTESTS_H_
#include "Fluxes.h"

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Select Numerical Scheme
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void (*calcFluxesPtr)();
void (*addRHSofRKPtr)(double ans[2], int j, int k);

/*
 * Calculate Fluxes: Wrapper Function
 */
void calcFluxes() {
	(*calcFluxesPtr)();
}

/*
 * Runge-Kutta Method RHS: Wrapper Function
 */
void addRHS_RK(double ans[2], int j, int k) {
	(*addRHSofRKPtr)(ans, j, k);
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Set Up The Test Parameters
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void setUpTests() {
	if (modelNumber == 0) {

	} else if (modelNumber == 1) {
		// Set up basics of test
		prep_Test1();
		// Set initial conditions
		initTFcnPtr = &exact_T_Test1;
		initqFcnPtr = &zeroInit;
		// Set source functions
		addSourceFcnPtr = &addSource_Test1;
		// Set boundary conditions
		dirichletBoundaryVal = 0;
		// Set flux functions
		if (numericalScheme == 0) {
			calcFluxesPtr = &calcFluxes_ClassFV_Dirichlet;
			addRHSofRKPtr = &addRHSofRK_FV;
		}
	}
}

#endif /* SETUPTESTS_H_ */
