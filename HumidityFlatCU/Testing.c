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
	prep_Test1();
	timeMethod();
	printf("\n");
	//printMatrix2DTimes2_Comp1(numCellsX, numCellsP, sl);
	// printMatrix2DTimes2_Comp1(numCellsX, numCellsP, Hx);
	printf("\nL2 norm of error = %f\n", relativeErrorL2norm());
	printf("\n");
}
