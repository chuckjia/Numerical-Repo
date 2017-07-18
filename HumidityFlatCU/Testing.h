/*
 * Testing.h
 *
 *  Created on: Jul 4, 2017
 *      Author: chuckjia
 */

#ifndef TESTING_H_
#define TESTING_H_
#include "TestFcns.h"


/* ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Testing
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

/* ----- ----- ----- ----- ----- -----
 * Mesh Testing
 * ----- ----- ----- ----- ----- ----- */

void printMesh_Test() {
	printf("Cell left coordinates\n");
	for (int i = 0; i < numCellsX; i++) {
		for (int j = 0; j < numCellsP; j++) {
			printf("%1.4f ", getCellLeftX(i, j));
		}
		printf("\n");
	}

	printf("Cell right coordinates\n");
	for (int i = 0; i < numCellsX; i++) {
		for (int j = 0; j < numCellsP; j++) {
			printf("%1.4f ", getCellRightX(i, j));
		}
		printf("\n");
	}

	printf("Cell bottom coordinates\n");
	for (int i = 0; i < numCellsX; i++) {
		for (int j = 0; j < numCellsP; j++) {
			printf("%1.4f ", getCellBottP(i, j));
		}
		printf("\n");
	}

	printf("Cell top coordinates\n");
	for (int i = 0; i < numCellsX; i++) {
		for (int j = 0; j < numCellsP; j++) {
			printf("%1.4f ", getCellTopP(i, j));
		}
		printf("\n");
	}
}

/*
 * Print u and omega cache for fluxes
 */
void printuomegaGodunov() {
	printf("Cached u values\n");
	for (int i = 0; i < numCellsX; i++) {
		for (int j = 0; j < numCellsP; j++) {
			printf("%1.4f ", uomega_Cache_Godunov[i][j][0] * DpInv);
		}
		printf("\n");
	}
	printf("Cached omega values\n");
	for (int i = 0; i < numCellsX; i++) {
		for (int j = numCellsP - 2; j < numCellsP; j++) {
			double x = getCellCenterX(i, j), p = getCellCenterP(i, j);
			double omegaExactVal = omega_fcn(x, p);
			printf("At (%1.5f, %1.5f), omega = %f; ", x, p, omegaExactVal);
			// printf("%1.10f ", uomega_Cache_Godunov[i][j][1] * DxInv);
		}
		printf("\n");
	}
}


/*
 * Peak value comparison for tests
 */

void compPeakValue_Test1() {
	// Print message
	printf("\nFunction value error:\n");
	double t = finalTime;

	// Location section units
	double numTestDivisions = 8;
	double xLocUnit = (double) Nx / numTestDivisions, pLocUnit = (double) Np / numTestDivisions;

	// Number of tests
	const int numCompSections = 4;
	double xLocFactor[numCompSections] = {1, 3, 5, 7},
			pLocFactor[numCompSections] = {1, 3, 5, 7};
	double maxError = 0, maxxLoc = 0, maxpLoc = 0;
	for (int ii = 0; ii < numCompSections; ii++)
		for (int jj = 0; jj < numCompSections; jj++) {
			int xLoc = (int) (xLocUnit * xLocFactor[ii]),
					pLoc = (int) (pLocUnit * pLocFactor[jj]);
			double x = getCellCenterX(xLoc, pLoc), p = getCellCenterP(xLoc, pLoc);
			double exactVal = (*initTFcnPtr)(x, p, t, 0, 0);
			double numericalVal = sl[xLoc][pLoc][0];
			double error = fabs(exactVal - numericalVal);
			printf(" - At section (%1.0f, %1.0f), coord (%1.1f, %1.0f): ",
					xLocFactor[ii], pLocFactor[jj], x * 1e-4, p);
			printf("error = %1.4f E-8, exact value = %f\n", error * 1e8, exactVal);
			if (error > maxError) {
				maxError = error;
				maxxLoc = xLocFactor[ii];
				maxpLoc = pLocFactor[jj];
			}
		}
	printf("Max error is = %1.4f E-8, Location (%1.0f, %1.0f)\n",
			maxError * 1e8, maxxLoc, maxpLoc);
}

#endif /* TESTING_H_ */
