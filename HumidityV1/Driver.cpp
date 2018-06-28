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
	setQuadCell();
	setProjection();
	setWPhix();
	setConditions();
	setAnalysis();
	setGodunov();
	setTimeSteps();
}

void fileManagement() {
	closeGlobalFiles_IO();
}

int main() {
	setAll();
	testing();
	fileManagement();
}
