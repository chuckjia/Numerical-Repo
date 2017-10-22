/*
 * Godunov.h
 *
 *  Created on: Oct 17, 2017
 *      Author: chuckjia
 */

#ifndef GODUNOV_H_
#define GODUNOV_H_
#include "Conditions.h"

// FF's are the fluxes on the right side of cell (i, j), where 0 <= i, j <= Nx (or lastRealX)
double FF_T[Nx + 1][Np + 1], FF_q[Nx + 1][Np + 1], FF_u[Nx + 1][Np + 1];
double GG_T[Nx + 1][Np + 1], GG_q[Nx + 1][Np + 1], GG_u[Nx + 1][Np + 1];

double r_uInterp_cache[Nx + 1][Np + 1];

void (*calcFluxes)();

void fillCache_r_uInterp_upwind() {
	for (int i = 0; i <= Nx; ++i)
		for (int j = 1; j <= Np; ++j) {
			double xMid = getCellRightX(i),
					leftPortion = xMid - getCellCenterX(i),
					rightPortion = getCellCenterX(i + 1) - xMid;
			r_uInterp_cache[i][j] = leftPortion / (leftPortion + rightPortion);
		}
}

double get_r_uInterp_upwind(int i, int j) {
	return r_uInterp_cache[i][j];
}

void upwind() {
	// GG fluxes: fluxes on TOP side of cells
	for (int i = 1; i <= Nx; ++i)
		for (int j = 0; j <= Np; ++j) {
			double uTopSideVal = 0.5 * (u[i][j] + u[i][j + 1]),
					wTopSideVal = 0.5 * (w[i][j] + w[i][j + 1]);
			double velocity = getCellTopSideLen() * (
					getCellTopSideNormX(i, j) * uTopSideVal +
					getCellTopSideNormP(i, j) * wTopSideVal);
			if (velocity >= 0) {
				GG_T[i][j] =	 velocity * T[i][j];
				GG_q[i][j] =	 velocity * q[i][j];
				GG_u[i][j] =	 velocity * u[i][j];
			} else {
				GG_T[i][j] =	 velocity * T[i][j + 1];
				GG_q[i][j] =	 velocity * q[i][j + 1];
				GG_u[i][j] =	 velocity * u[i][j + 1];
			}
		}

	// FF fluxes: fluxes on RIGHT side of cells
	for (int i = 0; i <= Nx; ++i)
		for (int j = 1; j <= Np; ++j) {
			double r = get_r_uInterp_upwind(i, j);
			double velocity = getCellRightSideLen(i, j) *
					(r * u[i + 1][j] + (1 - r) * u[i][j]);
			if (velocity >= 0) {
				GG_T[i][j] =	 velocity * T[i][j];
				GG_q[i][j] =	 velocity * q[i][j];
				GG_u[i][j] =	 velocity * u[i][j];
			} else {
				GG_T[i][j] =	 velocity * T[i][j + 1];
				GG_q[i][j] =	 velocity * q[i][j + 1];
				GG_u[i][j] =	 velocity * u[i][j + 1];
			}
		}
}

void setGodunov() {
	// if (fluxMethod == 0) {
	fillCache_r_uInterp_upwind();
	calcFluxes = &upwind;
}

#endif /* GODUNOV_H_ */
