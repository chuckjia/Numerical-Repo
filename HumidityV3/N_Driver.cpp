/*
 * Driver.cpp
 *
 *  Created on: Oct 13, 2017
 *      Author: chuckjia
 */

#include "M_Testing.h"

void setAll() {
	printTitle();
	setAllSettings();
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

void fileManagement() {
	closeFiles_analysis();
}

int main() {
	setAll();
	printCellCentersToFile();
	printParamToFile();
	runTimeSteps();
	showL2Errors();
}
