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

double relativeErrorL2norm() {
	double t = finalTime;
	double sum1 = 0;
	double sum2 = 0;
	for (int j = 1; j < lastIndexX; j++)
		for (int k = 1; k < lastIndexP; k++) {
			double x = getCellCenterX(j, k), p = getCellCenterP(j, k);
			double ExactVal = soln_T_Test1(x, p, t, j, k);
			double NumericalVal = sl[j][k][0];
			sum1 += pow(ExactVal - NumericalVal, 2);
			sum2 += ExactVal * ExactVal;
		}
	return sqrt(sum1 / sum2);
}

double absErrorL2norm() {
	double t = finalTime;
	double sum1 = 0;
	double sum2 = 0;
	for (int j = 1; j < lastIndexX; j++)
		for (int k = 1; k < lastIndexP; k++) {
			double x = getCellCenterX(j, k), p = getCellCenterP(j, k);
			double ExactVal = soln_T_Test1(x, p, t, j, k);
			double NumericalVal = sl[j][k][0];
			sum1 += pow(ExactVal - NumericalVal, 2);
		}
	return sqrt(cellVol * sum1);
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
