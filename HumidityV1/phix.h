/*
 * phix.h
 *
 *  Created on: Oct 23, 2017
 *      Author: chuckjia
 */

#ifndef PHIX_H_
#define PHIX_H_
#include "w.h"

double getCellTopRightT(int i, int j) {
	return a1_quadCell * u[i][j] + a2_quadCell[i][j] * u[i + 1][j] +
			a3_quadCell[i][j] * u[i][j + 1] + a4_quadCell[i][j] * u[i + 1][j + 1];
}

double getTGradX(int i, int j) {
	return e11_diagMatInv_quadCell * (getCellTopRightT(i, j) - getCellTopRightT(i - 1, j)) +
			e12_diagMatInv_quadCell[i][j] * (T[i][j + 1] - T[i][j]);
}

double getTGradP(int i, int j) {
	return e22_diagMatInv_quadCell[i][j] * (T[i][j + 1] - T[i][j]);
}



#endif /* PHIX_H_ */
