/*
 * WPhix.h
 *
 *  Created on: Oct 21, 2017
 *      Author: chuckjia
 *
 *  Functions that computes w and phi_x.
 */

#ifndef G_WPHIX_H_
#define G_WPHIX_H_
#include "F_Projection.h"

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Calculate the w Function
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void (*calcW_fptr)();

// Calculate values of w
void calcW_orig() {
	for (int i = 1; i <= Nx; ++i) {
		// w_[i][0] = -w_[i][1];  // Might be unnecessary, just to be safe
		for (int j = 0; j < Np; ++j)
			w_[i][j + 1] = w_[i][j] - (getCellCenterP(i, j + 1) - getCellCenterP(i, j)) * getGradhU_x(i, j);
	}
}


/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Calculate the phi_x Function
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void (*calcPhix_fptr)();

// This function uses formula in (3.42): questionable calculation
void calcPhix_orig() {
	for (int i = 1; i <= Nx; ++i) {
		//		phix_[i][0] = 0;  // Might be unnecessary, just to be safe
		//		phix_[i][1] = 0;
		double sum = 0, DpVal = getCellCenterDp(i);
		for (int j = 1; j < Np; ++j) {
			// sum += (getCellTopCenterP(i, j) - getCellBottCenterP(i, j)) / getCellCenterP(i, j) * getGradhT_x(i, j);  // Changed! Back to old implementation
			sum += DpVal * getGradhT_x(i, j) / getCellCenterP(i, j);
			phix_[i][j + 1] = -R_CONST * sum;
		}
	}
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Set All Parameters
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void setWPhix() {
	calcW_fptr = &calcW_orig;
	calcPhix_fptr = &calcPhix_orig;
}

#endif /* G_WPHIX_H_ */
