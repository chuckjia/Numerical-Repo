/*
 * w.h
 *
 *  Created on: Oct 21, 2017
 *      Author: chuckjia
 */

#ifndef WPHIX_H_
#define WPHIX_H_
#include "QuadCell.h"

void calc_w() {
	for (int i = 1; i <= Nx; ++i)
		for (int j = 0; j < Np; ++j)
			w[i][j + 1] = w[i][j] - (getCellCenterP(i, j + 1) - getCellCenterP(i, j)) *
					getGradhU_x(i, j);
}

void calc_phix() {
	for (int i = 1; i <= Nx; ++i) {
		double sum = 0;
		// factor = -R * (p_{i,j+1/2} - p_{i,j-1/2})
		double factor = -0.5 * R_CONST * (getCellLeftDp(i) + getCellRightDp(i));
		for (int j = 1; j < Np; ++j) {
			sum += factor / getCellCenterP(i, j) * getGradhT_x(i, j);
			phix[i][j] = sum;
		}
	}
}

#endif /* WPHIX_H_ */
