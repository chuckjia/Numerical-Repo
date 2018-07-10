/*
 * Godunov.h
 *
 *  Created on: Oct 17, 2017
 *      Author: chuckjia
 *
 *  Methods and arrays to apply the Godunov upwind method
 */

#ifndef J_GODUNOV_H_
#define J_GODUNOV_H_
#include "I_IO.h"

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Interpolation for u
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Cache to store interpolation coefficient r used in (3.31) for u_{i+1/2,j} values. r[i][j] is used to interpolate for the u value on the
// right side of cell(i,j). The value of r is invariant with respect to j, by using the similarity property of two triangles
double r_uInterp_[Nx + 1];

// Calculate r values
void fillCache_r_uInterp() {
	for (int i = 0; i <= Nx; ++i){
		// The naming corresponds to Figure 6 on page 110
		double xCommonBD = getCellRightX(i),
				leftPortion = xCommonBD - getCellCenterX(i),
				rightPortion = getCellCenterX(i + 1) - xCommonBD;
		r_uInterp_[i] = leftPortion / (leftPortion + rightPortion);
	}
}

// Return r value for the calculation of u_{i+1/2,j}
double get_r_uInterp(int i, int j) { return r_uInterp_[i]; }
double get_r_uInterp(int i) { return r_uInterp_[i]; }


/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Fluxes and Flux Calculation
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// FF: fluxes on the RIGHT side of cell (i, j), with 0 <= i <= Nx, 0 <= j <= Np
double FF_T_[Nx + 1][Np + 1], FF_q_[Nx + 1][Np + 1], FF_u_[Nx + 1][Np + 1];  // FF_{i+1/2,j}
// GG: fluxes on the TOP side of cell (i, j), with 0 <= i <= Nx, 0 <= j <= Np
double GG_T_[Nx + 1][Np + 1], GG_q_[Nx + 1][Np + 1], GG_u_[Nx + 1][Np + 1];  // GG_{i,j+1/2}

// Function pointer: function to calculate all fluxes
void (*calcFluxes)();

// Calculate all fluxes by upwind type Godunov scheme
void calcFluxes_upwind() {
	// GG: fluxes on the TOP sides of cells
	for (int i = 1; i <= Nx; ++i)
		for (int j = 0; j <= Np; ++j) {
			double uTopSideVal = 0.5 * (u_[i][j] + u_[i][j + 1]),
					wTopSideVal = 0.5 * (w_[i][j] + w_[i][j + 1]);
			double velocity = getCellTopSideLen(i, j) * (
					getCellTopSideNormVecX(i, j) * uTopSideVal + getCellTopSideNormVecP(i, j) * wTopSideVal);
			if (velocity >= 0) {
				GG_T_[i][j] = velocity * T_[i][j];
				GG_q_[i][j] = velocity * q_[i][j];
				GG_u_[i][j] = velocity * u_[i][j];
			} else {
				GG_T_[i][j] = velocity * T_[i][j + 1];
				GG_q_[i][j] = velocity * q_[i][j + 1];
				GG_u_[i][j] = velocity * u_[i][j + 1];
			}
		}

	// FF: fluxes on the RIGHT sides of cells
	for (int i = 0; i <= Nx; ++i) {
		double r = get_r_uInterp(i);
		for (int j = 1; j <= Np; ++j) {
			double velocity = getCellRightSideLen(i, j) * (r * u_[i + 1][j] + (1 - r) * u_[i][j]);
			if (velocity >= 0) {
				FF_T_[i][j] = velocity * T_[i][j];
				FF_q_[i][j] = velocity * q_[i][j];
				FF_u_[i][j] = velocity * u_[i][j];
			} else {
				FF_T_[i][j] = velocity * T_[i + 1][j];
				FF_q_[i][j] = velocity * q_[i + 1][j];
				FF_u_[i][j] = velocity * u_[i + 1][j];
			}
		}
	}
}


/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Set All Parameters for Godunov Type Flux Calculation
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void setGodunov() {
	fillCache_r_uInterp();
	calcFluxes = &calcFluxes_upwind;
}

#endif /* J_GODUNOV_H_ */
