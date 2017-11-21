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
	setQuadCells();
	setProjection();
	setWPhix();
	setConditions();
	setAnalysis();
	setGodunov();
	setTimeSteps();
}

int main() {
	setAll();
	testing();
}
