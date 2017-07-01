/*
 * Driver.c
 *
 *  Created on: Jun 10, 2017
 *      Author: chuckjia
 */

#include "TimeMethod.h"

int main() {
	buildMesh();
	setInitCond();
	//printMatrix3D_1(numCellsXDir, numCellsPDir, soln);
	rk2_CU();
	//printMatrix3D_1(numCellsXDir, numCellsPDir, soln);
}
