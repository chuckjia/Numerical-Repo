/*
 * Driver.cpp
 *
 *  Created on: Jun 20, 2017
 *      Author: chuckjia
 */

#include "Testing.h"

void math();
void test();

int main() {
	math();
	test();
}

void math() {
	printf("\n===== ===== ===== ===== ===== ===== ===== ===== ===== ===== \n");
	printf(": Solving Model Equations\n");
	printf("===== ===== ===== ===== ===== ===== ===== ===== ===== ===== \n");

	printf("\nProgram Progress\n");
	setMesh();
	setUpTests();
	setInitCond();
	timeMethod();
	// printMatrix2DTimes2_Comp1(numCellsX, numCellsP, sl);
	// printMatrix2DTimes2_Comp1(numCellsX, numCellsP, Hx);
	printf("\n----- ----- ----- ----- ----- ----- ----- -----\n");
	printf("\nMethod: Classical FV");
	printf("\nSpace:  Nx = %d | Np = %d || Dx = %1.2f, Dp = %1.2f\n",
			Nx, Np, Dx, Dp);
	printf("\nTime:   Steps = %d | Final Time = %1.4fs\n", numTimeSteps, finalTime);
	printf("\nModel:  Test %d\n", modelNumber);
	if (modelNumber != 0) {
		printf("\n----- ----- ----- ----- ----- ----- ----- -----\n");
		printf("\nComparison Data\n");
		showL2Errors();
		showL1Errors();
		middleDiff();
	}
	writeResToFile();
	printf("\n");
}

void test() {
}
