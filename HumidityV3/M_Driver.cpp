/*
 * Driver.cpp
 *
 *  Created on: Oct 13, 2017
 *      Author: chuckjia
 */

#include "L_Testing.h"

void setAll() {
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

int main() {
	setAll();
	writeCSV_cellCenters();
	writeCSV_param();
	runTimeSteps();
	showL2Errors();
}
