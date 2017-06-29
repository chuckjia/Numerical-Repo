/*
 * Driver.c
 *
 *  Created on: Jun 20, 2017
 *      Author: chuckjia
 */

#include "TestFcns.h"

int main() {
	printf("\n===== ===== ===== ===== ===== ===== ===== ===== ===== ===== \n");
	printf("Solving model equations using ");
	if (numericalScheme == 1)
		printf("the Central Upwind Method\n");
	else
		printf("the Godnuv method\n");
	printf("===== ===== ===== ===== ===== ===== ===== ===== ===== ===== \n");
	printf("Parameters:\n  Nx = %d, Np = %d, numTimeSteps = %d\n", Nx, Np, numTimeSteps);

	buildMesh();
	setUpTests();
	setInitCond();
	selectScheme();
	timeMethod();
	// printMatrix2DTimes2_Comp1(numCellsX, numCellsP, sl);
	// printMatrix2DTimes2_Comp1(numCellsX, numCellsP, Hx);
	printf("\nL2 norm of relative error = %1.10f\n", relativeErrorL2norm());
	printf("\nL2 norm of absolute error = %1.10f\n", absErrorL2norm());
	writeResToFile();
	printf("\n");
}
