/*
 * w.h
 *
 *  Created on: Oct 21, 2017
 *      Author: chuckjia
 */

#ifndef WPHIX_H_
#define WPHIX_H_
#include "QuadCell.h"

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Calculate the w Function
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void calc_w() {
	for (int i = 1; i <= Nx; ++i)
		for (int j = 0; j < Np; ++j)
			w_sl[i][j + 1] = w_sl[i][j] - (getCellCenterP(i, j + 1) - getCellCenterP(i, j)) *
					getGradhU_x(i, j);
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Calculate the phi_x Function
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// This function uses formula in (3.42): questionable calculation
void calc_phix() {
	for (int i = 1; i <= Nx; ++i) {
		double sum = 0;
		// factor = -R * (p_{i,j+1/2} - p_{i,j-1/2})
		double factor = -R_CONST * getCellCenterDp(i);
		for (int j = 1; j < Np; ++j) {
			sum += factor / getCellCenterP(i, j) * getGradhT_x(i, j);
			phix_sl[i][j + 1] = sum;
		}
	}
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Set All Parameters
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void setWPhix() {
	// Empty for now
}

#endif /* WPHIX_H_ */
