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

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Testing
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void testing() {
	enforceIC();
	print_lambdax();
}
#endif /* L_TESTING_H_ */
