/*
 * TestFcns.h
 *
 *  Created on: Jul 17, 2017
 *      Author: chuckjia
 */

#ifndef TESTFCNS_H_
#define TESTFCNS_H_
#include "TimeSteps.h"

void showL2Errors() {
	double t = finalTime;
	double sum1 = 0, sum2 = 0;
	for (int i = 1; i < lastGhostIndexX; i++)
		for (int j = 1; j < lastGhostIndexP; j++) {
			double x = getCellCenterX(i, j), p = getCellCenterP(i, j);
			double ExactVal = exact_T_Test1(x, p, t);
			double NumericalVal = soln[i][j];
			sum1 += pow(ExactVal - NumericalVal, 2);
			sum2 += ExactVal * ExactVal;
		}
	double absError = sqrt(cellVol * sum1);
	double relativeError = 0;
	if (sum1 > 1e-15)
		relativeError = sqrt(sum1 / sum2);
	printf("\n- L2 relative error = %1.15f\n", relativeError);
	printf("\n- L2 absolute error = %1.15f\n", absError);
}

void printMsg() {
	printf("\nnumDivisions = %d, numTimeSteps = %d, Dt = %f, finalTime = %1.4f\n",
			numDivisions, numTimeSteps, Dt, finalTime);
}

void writeResToFile() {
	FILE *f = fopen("res.txt", "wb");
	for (int i = 1; i < lastGhostIndexX; i++)
		for (int j = 1; j < lastGhostIndexP; j++)
			fprintf(f, "%f ", soln[i][j]);
	fclose(f);
	FILE *g = fopen("err.txt", "wb");
	for (int i = 1; i < lastGhostIndexX; i++)
		for (int j = 1; j < lastGhostIndexP; j++) {
			double x = getCellCenterX(i, j), p = getCellCenterP(i, j);
			double exactSoln = exact_T_Test1(x, p, finalTime);
			fprintf(g, "%f ", soln[i][j] - exactSoln);
		}
	fclose(g);
}

void printMat(int numRows, int numCols, double mat[numRows][numCols]) {
	for (int i = 0; i < numRows; i++) {
		for (int j = 0; j < numCols; j++) {
			printf("%1.2f ", mat[i][j]);
		}
		printf("\n");
	}
}

void printFlux(double mat[numCellsX][numCellsP]) {
	for (int i = 1; i < lastGhostIndexX; i++) {
		for (int j = 1; j < lastGhostIndexP; j++) {
			printf("%1.2f ", mat[i][j]);
		}
		printf("\n");
	}
}

#endif /* TESTFCNS_H_ */
