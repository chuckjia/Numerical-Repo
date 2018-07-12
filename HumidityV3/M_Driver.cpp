/*
 * Driver.cpp
 *
 *  Created on: Oct 13, 2017
 *      Author: chuckjia
 */

#include "L_Testing.h"

void setAll() {
	if (movieFrameFreq <= 0)
		movieFrameFreq = numTimeStep + 1;
	if (calcL2errFreq <= 0)
		calcL2errFreq = numTimeStep + 1;

	printTitle();
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
	writeCSV_param();
	runTimeSteps();
	showL2Errors();
	writeCSV_finalSolnErr();
	printf("\n");
}
