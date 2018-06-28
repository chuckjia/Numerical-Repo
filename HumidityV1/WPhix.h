/*
 * w.h
 *
 *  Created on: Oct 21, 2017
 *      Author: chuckjia
 */

#ifndef WPHIX_H_
#define WPHIX_H_
#include "Projection.h"

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Calculate the w Function
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void (*calcW_fptr)();

void calcW_orig() {
	for (int i = 1; i <= Nx; ++i)
		for (int j = 0; j < Np; ++j)
			w_[i][j + 1] = w_[i][j] - (getCellCenterP(i, j + 1) - getCellCenterP(i, j)) *
					getGradhU_x(i, j);
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Calculate the phi_x Function
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void (*calcPhix_fptr)();

// This function uses formula in (3.42): questionable calculation
void calcPhix_orig() {
	for (int i = 1; i <= Nx; ++i) {
		double sum = 0;
		// factor = -R * (p_{i,j+1/2} - p_{i,j-1/2})
		double factor = -R_CONST * getCellCenterDp(i);
		for (int j = 1; j < Np; ++j) {
			sum += factor / getCellCenterP(i, j) * getGradhT_x(i, j);
			phix_[i][j + 1] = sum;
		}
		//phix_sl[i][1] = phix_sl[i][2];
	}
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Set All Parameters
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void setWPhix() {
	calcW_fptr = &calcW_orig;
	calcPhix_fptr = &calcPhix_orig;
}

#endif /* WPHIX_H_ */
