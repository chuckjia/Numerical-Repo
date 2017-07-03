/*
 * Driver.c
 *
 *  Created on: Jun 20, 2017
 *      Author: chuckjia
 */

#include "TestFcns.h"

int main() {
	printf("\n===== ===== ===== ===== ===== ===== ===== ===== ===== ===== \n");
	if (numericalScheme == 1)
		printf("Central Upwind Method");
	else
		printf("Godunov Method");
	printf(": Solving Model Equations\n");
	printf("===== ===== ===== ===== ===== ===== ===== ===== ===== ===== \n");
	printf("\nParameters:  Nx = %d, Np = %d, numTimeSteps = %d, finalTime = %1.2fs\n",
			Nx, Np, numTimeSteps, finalTime);
	printf("Model Selection: Test %d\n", testNumber);
	printf("\n----- ----- ----- -----\n");
	printf("\nProgram Progress\n");
	buildMesh();
	setUpTests();
	setInitCond();
	selectScheme();
	timeMethod();
	// printMatrix2DTimes2_Comp1(numCellsX, numCellsP, sl);
	// printMatrix2DTimes2_Comp1(numCellsX, numCellsP, Hx);
	if (testNumber != 0) {
		printf("\n----- ----- ----- -----\n");
		printf("\nComparison Data\n");
		showL2Errors();
		showL1Errors();
		middleDiff();
	}
	writeResToFile();
	printf("\n");
}
