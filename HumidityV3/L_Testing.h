/*
 * Testing.h
 *
 *  Created on: Oct 16, 2017
 *      Author: chuckjia
 */

#ifndef L_TESTING_H_
#define L_TESTING_H_
#include "K_TimeSteps.h"

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Testing the LU Solver
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void projU_diagnostics() {
	/*printf("\n");
	double prevSum = 0;
	for (int i = 1; i <= Nx; ++i)
		if (!(i % 10)) {
			double currSum = 0;
			for (int j = 1; j <= Np; ++j)
				currSum += u_[i][j];
			currSum *= getCellCenterDp(i);
			// printf("Int u[%d] = %1.4f,  ", i, currSum);
			printf("Dx Int u[%d] = %1.4f,  ", i, (currSum - prevSum) / Dx);
			prevSum = currSum;
		}
	printf("\n");*/

	for (int i = 0; i <= Nx; ++i) {
		double x = getCellCenterX(i);
		double pBxDer = (*pBxDer_fptr)(x);
		printf("i=%d, w(pB) - u(pB)*(DpB/Dx) = %e,    ", i, w_[i][Np] - u_[i][Np] * pBxDer);
	}
	printf("\n");
}

//void test_baryCenter(int i, int j) {
//	double center[2];
//	double x1 = getCellLeftX(i, j), x2 = getCellRightX(i, j), x3 = x2, x4 = x1;
//	double p1 = getCellBottLeftP(i, j), p2 = getCellBottRightP(i, j), p3 = getCellTopRightP(i, j), p4 = getCellTopLeftP(i, j);
//	calcTrapezoidCenter(center, x1, x2, x3, x4, p1, p2, p3, p4);
//	printf(">> x_left = %1.5e, x_right = %1.5e\n", x1, x2);
//	printf(">> p_bottLeft = %1.5e, p_bottRight = %1.5e, p_topRight = %1.5e, p_topLeft = %1.5e\n", p1, p2, p3, p4);
//	printf(">> x_center =                   p_center = \n");
//	printf("   %1.20e,  %1.20e\n\n", center[0], center[1]);
//}

void writeCSV_testQuad() {
	FILE *fu = fopen("Output/gradhU_x.csv", "wb");
	int last_j = Np - 1;
	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < last_j; ++j)
			fprintf(fu, "%1.20e,", getGradhU_x(i, j));
		fprintf(fu, "%1.20e\n", getGradhU_x(i, last_j));
	}
	fclose(fu);
}

void writeCSV_testQuadCoef_helper(double mat[Nx + 1][Np + 1], string filename) {
	FILE *f = fopen(filename.c_str(), "wb");
	for (int i = 0; i <= Nx; ++i) {
		for (int j = 0; j < Np; ++j)
			fprintf(f, "%1.20e,", mat[i][j]);
		fprintf(f, "%1.20e\n", mat[i][Np]);
	}
	fclose(f);
}

void writeCSV_testQuadCoef() {
	printf(">> Printed all quad coefficients a2, a3, and a4 to CSV file.\n");
	writeCSV_testQuadCoef_helper(a2_quad_, "Output/a2_quad.csv");
	writeCSV_testQuadCoef_helper(a3_quad_, "Output/a3_quad.csv");
	writeCSV_testQuadCoef_helper(a4_quad_, "Output/a4_quad.csv");
	writeCSV_testQuadCoef_helper(e12_MInv_quad_, "Output/e12_MInv_quad.csv");

}

void writeCSV_testCellTopRightU() {
	FILE *f = fopen("Output/CellTopRightU.csv", "wb");
	for (int i = 0; i <= Nx; ++i) {
		for (int j = 0; j < Np; ++j)
			fprintf(f, "%1.20e,", getCellTopRightU(i, j));
		fprintf(f, "%1.20e\n", getCellTopRightU(i, Np));
	}
	fclose(f);
}

void testSource() {
	double x = 10.3847, p = 10.567, t = 20.123;
	printf("At x = %1.5e, p = %1.5e, t = %1.5e,\n  Source = %1.25e\n", x, p, t, (*source_T_fptr)(1, 1, 1, 1, x, p, t));
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Testing
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void testing() {
	runTimeSteps();
}

#endif /* L_TESTING_H_ */
