/*
 * TestFcns.h
 *
 *  Created on: Jun 20, 2017
 *      Author: chuckjia
 */

#ifndef TESTFCNS_H_
#define TESTFCNS_H_
#include "TimeSteps.h"

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Methods For Testing Purposes
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void showL2Errors() {
	double t = finalTime;
	double sum1 = 0;
	double sum2 = 0;
	for (int j = 1; j < lastGhostX; j++)
		for (int k = 1; k < lastGhostP; k++) {
			double x = getCellCenterX(j, k), p = getCellCenterP(j, k);
			double ExactVal = (*initTFcnPtr)(x, p, t, j, k);
			double NumericalVal = soln[j][k][0];
			sum1 += pow(ExactVal - NumericalVal, 2);
			sum2 += ExactVal * ExactVal;
		}
	double absError = sqrt(cellVol * sum1);
	double relativeError = sqrt(sum1 / sum2);
	if (absError < 1e-16)
		relativeError = 0;
	printf("\n- L2 relative error = %1.15f\n", relativeError);
	printf("\n- L2 absolute error = %1.15f\n", absError);
}

void showL1Errors() {
	double t = finalTime;
	double sum1 = 0;
	double sum2 = 0;
	for (int j = 1; j < lastGhostX; j++)
		for (int k = 1; k < lastGhostP; k++) {
			double x = getCellCenterX(j, k), p = getCellCenterP(j, k);
			double ExactVal = (*initTFcnPtr)(x, p, t, j, k);
			double NumericalVal = soln[j][k][0];
			sum1 += fabs(ExactVal - NumericalVal);
			sum2 += fabs(ExactVal);
		}
	printf("\n- L1 relative error = %1.8f E-3\n", sqrt(sum1 / sum2) * 1e3);
	printf("\n- L1 absolute error = %1.15f\n", sqrt(cellVol * sum1));
}

void middleDiff() {
	double xMid = x0 + 0.5 * (xf - x0), pMid = pA + 0.5 * (pB - pA);
	int j = Nx / 2, k = Np / 2;
	double exactSolnCenter = (*initTFcnPtr)(xMid, pMid, finalTime, j, k),
			numericalSolnCenter = soln[j][k][0];
	printf("\n- Center value comparison\n");
	printf("    Error = %f\n", fabs(exactSolnCenter - numericalSolnCenter));
	printf("    Exact value = %f", exactSolnCenter);
}

void writeResToFile() {
	FILE *f = fopen("resT.txt", "wb");
	for (int i = 0; i < numCellsX; i++) {
		for (int j = 0; j < numCellsP; j++) {
			double val = soln[i][j][0];
			if (fabs(val) < 0.00001)
				val = 0;
			fprintf(f, "%f ", val);
		}
		fprintf(f, "\n");
	}
	fclose(f);
	FILE *g = fopen("resq.txt", "wb");
	for (int i = 0; i < numCellsX; i++) {
		for (int j = 0; j < numCellsP; j++) {
			double val = soln[i][j][1];
			if (fabs(val) < 0.00001)
				val = 0;
			fprintf(g, "%f ", val);
		}
		fprintf(g, "\n");
	}
	fclose(g);
}

#endif /* TESTFCNS_H_ */
