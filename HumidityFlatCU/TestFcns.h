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
	for (int j = 1; j < lastIndexX; j++)
		for (int k = 1; k < lastIndexP; k++) {
			double x = getCellCenterX(j, k), p = getCellCenterP(j, k);
			double ExactVal = (*calcExactT)(x, p, t, j, k);
			double NumericalVal = sl[j][k][0];
			sum1 += pow(ExactVal - NumericalVal, 2);
			sum2 += ExactVal * ExactVal;
		}
	printf("\n- L2 relative error = %1.10f\n", sqrt(sum1 / sum2));
	printf("\n- L2 absolute error = %1.10f\n", sqrt(cellVol * sum1));
}

void showL1Errors() {
	double t = finalTime;
	double sum1 = 0;
	double sum2 = 0;
	for (int j = 1; j < lastIndexX; j++)
		for (int k = 1; k < lastIndexP; k++) {
			double x = getCellCenterX(j, k), p = getCellCenterP(j, k);
			double ExactVal = (*calcExactT)(x, p, t, j, k);
			double NumericalVal = sl[j][k][0];
			sum1 += fabs(ExactVal - NumericalVal);
			sum2 += fabs(ExactVal);
		}
	printf("\n- L1 relative error = %1.10f\n", sqrt(sum1 / sum2));
	printf("\n- L1 absolute error = %1.10f\n", sqrt(cellVol * sum1));
}

void middleDiff() {
	double xMid = x0 + 0.5 * (xL - x0), pMid = pA + 0.5 * (pB - pA);
	int j = Nx / 2, k = Np / 2;
	double exactSolnCenter = (*calcExactT)(xMid, pMid, finalTime, j, k),
			numericalSolnCenter = sl[j][k][0];
	printf("\n- Center value (peak) comparison\n");
	printf("    Error = %f\n", fabs(exactSolnCenter - numericalSolnCenter));
	printf("    Exact value = %f\n", exactSolnCenter);

}

void writeResToFile() {
	FILE *f = fopen("res.txt", "wb");
	for (int i = 0; i < numCellsX; i++) {
		for (int j = 0; j < numCellsP; j++) {
			double val = sl[i][j][0];
			if (fabs(val) < 0.00001)
				val = 0;
			fprintf(f, "%f ", val);
		}
		fprintf(f, "\n");
	}
	fclose(f);
}

#endif /* TESTFCNS_H_ */
