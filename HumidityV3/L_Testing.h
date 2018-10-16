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

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Testing
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void testing() {
	//	writeCSV_quadCoefs();
	//	writeCSV_matrix(w_, "Output/w_test.csv");
	//	writeCSV_MInv_quad();
	//	test_u_quad();
	//	test_gradhUx();
	runTimeSteps();
}

#endif /* L_TESTING_H_ */
