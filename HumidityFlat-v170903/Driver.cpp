/*
 * Driver.c
 *
 *  Created on: Jul 19, 2017
 *      Author: chuckjia
 */

#include "FV.h"

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Result Analysis
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Print the L2 relative and absolute errors
void showL2Errors() {
	double t = finalTime;
	double sum1 = 0, sum2 = 0;
	for (int i = 1; i < lastGhostIndexX; i++)
		for (int j = 1; j < lastGhostIndexP; j++) {
			double x = getCellCenterX(i, j), p = getCellCenterP(i, j);
			double exactVal = (*initTPtr)(x, p, t);
			double numericalVal = sl[i][j][0];
			sum1 += pow(exactVal - numericalVal, 2);
			sum2 += exactVal * exactVal;
		}
	double absError = sqrt(cellVol * sum1);
	double relativeError = 0;
	if (sum1 > 1e-15)
		relativeError = sqrt(sum1 / sum2);
	printf("\n- L2 relative error = %1.15f\n", relativeError);
	printf("\n- L2 absolute error = %1.15f\n", absError);
}

// Print the L2 relative and absolute errors
void showL1Errors() {
	double t = finalTime;
	double sum1 = 0, sum2 = 0;
	for (int i = 1; i < lastGhostIndexX; i++)
		for (int j = 1; j < lastGhostIndexP; j++) {
			double x = getCellCenterX(i, j), p = getCellCenterP(i, j);
			double exactVal = (*initTPtr)(x, p, t);
			double numericalVal = sl[i][j][0];
			sum1 += fabs(exactVal - numericalVal);
			sum2 += fabs(exactVal);
		}
	double absError = cellVol * sum1;
	double relativeError = 0;
	if (sum1 > 1e-15)
		relativeError = sum1 / sum2;
	printf("\n- L1 relative error = %1.15f\n", relativeError);
	printf("\n- L1 absolute error = %1.15f\n", absError);
}

// Check which boundary conditions are used
void printBC() {
	if (enforceBCPtr == &enforcePeriodicBC)
		printf(", BC: Periodic");
	else if (enforceBCPtr == &enforceDirichletBC)
		printf(", BC: Dirichlet");
	else if (enforceBCPtr == &enforceNeumannBC)
		printf(", BC: Neumann");
}

// Print out messages of scheme information in the console
void printMsg() {
	printf("\nModel Number: %d, ", modelNum);
	if (fluxMethod == 0)
		printf("Method: Godunov");
	else
		printf("Method: Classical FV");
	printBC(); printf("\n");
	printf("\nnumTimeSteps = %d, Dt = %1.2e, finalTime = %1.4f\n", numTimeSteps, Dt, finalTime);
	printf("numDivisions = %d, Dx = %1.2e, Dt/Dx = %1.2f\n", numDivisions, Dx, Dt / Dx);
}

// Write the result to file: res.txt for the solution and err.txt for the error
// This function is used in ploting the solution and the error
void writeResToFile() {
	FILE *f = fopen("resT.txt", "wb");
	for (int i = 1; i < lastGhostIndexX; i++)
		for (int j = 1; j < lastGhostIndexP; j++)
			fprintf(f, "%f ", sl[i][j][0]);
	fclose(f);
	FILE *g = fopen("errT.txt", "wb");
	for (int i = 1; i < lastGhostIndexX; i++)
		for (int j = 1; j < lastGhostIndexP; j++) {
			double x = getCellCenterX(i, j), p = getCellCenterP(i, j);
			double exactSoln = (*initTPtr)(x, p, finalTime);
			fprintf(g, "%f ", sl[i][j][0] - exactSoln);
		}
	fclose(g);
}

void printErr() {
	for (int i = 1; i < lastGhostIndexX; i++) {
		for (int j = 1; j < lastGhostIndexP; j++){
			double x = getCellCenterX(i, j), p = getCellCenterP(i, j);
			double exactSoln = (*initTPtr)(x, p, finalTime);
			printf("%1.4f ", sl[i][j][0] - exactSoln);
		}
		printf("\n");
	}
}

void printAnalysis() {
	printMsg();
	showL2Errors();
	showL1Errors();
	writeResToFile();
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Main Function
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

int main() {
	selectModel();
	forwardEuler();
	printAnalysis();
}
