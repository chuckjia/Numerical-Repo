/*
 * Driver.c
 *
 *  Created on: Jun 20, 2017
 *      Author: chuckjia
 */

#include "TestFcns.h"

int main() {
	printf("===== ===== ===== ===== ===== ===== ===== ===== ===== ===== \n");
	printf("Godunov Method: Solving Model Equations\n");
	printf("===== ===== ===== ===== ===== ===== ===== ===== ===== ===== \n");
	printf("\nProgram Progress\n");
	buildMesh();
	setUpTests();
	setInitCond();
	timeMethod();
	// printMatrix2DTimes2_Comp1(numCellsX, numCellsP, sl);
	// printMatrix2DTimes2_Comp1(numCellsX, numCellsP, Hx);
	printf("\n----- ----- ----- -----\n");
	printf("\nParameters:  Nx = %d, Np = %d, numTimeSteps = %d, finalTime = %1.2fs\n",
			Nx, Np, numTimeSteps, finalTime);
	printf("Model Selection: Test %d\n", testNumber);
	printf("\n----- ----- ----- -----\n");
	printf("\nComparison Data\n");
	showL2Errors();
	showL1Errors();
	middleDiff();
	writeResToFile();
	printf("\n");
}
