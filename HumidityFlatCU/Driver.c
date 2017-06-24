/*
 * Driver.c
 *
 *  Created on: Jun 20, 2017
 *      Author: chuckjia
 */

#include "TimeSteps.h"

void normL2() {

}

int main() {
	buildMesh();
	setInitCond();
	rk2();
	printMatrix2DTimes2_Comp1(numCellsX, numCellsP, sl);
}
