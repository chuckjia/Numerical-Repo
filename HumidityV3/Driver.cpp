/*
 * Driver.cpp
 *
 *  Created on: Oct 13, 2017
 *      Author: chuckjia
 */

#include "L_Testing.h"

void setAll() {
	printTitle();
	validateProgramParameters();
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

int main(int argc, char *argv[]) {
	if (argc > 1)
		runInEclipse = false;
	setAll();
	printSchemeSummary();
	writeCSV_cellCenters();
	writeCSV_meshGrid();
	writeCSV_param();
	runTimeSteps();
	// testing();
	showL2Errors();
	writeCSV_finalSolnErr();
	printf("\n");
}
