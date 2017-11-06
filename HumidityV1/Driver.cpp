/*
 * Driver.cpp
 *
 *  Created on: Oct 13, 2017
 *      Author: chuckjia
 */

#include "Testing.h"

void setAll() {
	printTitle();
	setConstants();
	setModels();
	setMesh();
	setConditions();
	setGodunov();
	setProjection();
	setQuadCells();
	setWPhix();
	setTimeSteps();
}

int main() {
	setAll();
	testing();
}
