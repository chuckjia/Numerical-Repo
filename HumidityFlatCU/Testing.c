/*
 * Testing.c
 *
 *  Created on: Jun 20, 2017
 *      Author: chuckjia
 */

#include "TestFcns.h"

int main() {
	printf("\n===== ===== ===== ===== ===== ===== ===== ===== ===== \n");
	printf("Solving model equations by the central upwind method.\n");
	printf("===== ===== ===== ===== ===== ===== ===== ===== ===== \n");
	buildMesh();
	setInitCond();
	test1Prep();
	timeMethod();
	printf("\n");
	//printMatrix2DTimes2_Comp1(numCellsX, numCellsP, sl);
	// printMatrix2DTimes2_Comp1(numCellsX, numCellsP, Hx);
	printf("\nL2 norm of error = %f\n", errorL2norm());
	printf("\n");
}
