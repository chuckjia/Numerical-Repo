/*
 * FV.h
 *
 *  Created on: Aug 24, 2017
 *      Author: chuckjia
 */

#ifndef FV_H_
#define FV_H_
#include "Models.h"

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Fluxes
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// GG_Bott and GG_Top are fluxes along the vertical direction (p direction)
double GG[numCellsX][numCellsP][2];  // GG_{i,j+1/2}
// FF_Left and FF_Right are fluxes along the horizontal direction (x direction)
double FF[numCellsX][numCellsP][2];  // FF_{i+1/2,j}

// Calculate numerical fluxes using Godunov method (upwind type)
void calcFluxes_Godunov() {
	// GG fluxes: fluxes on top sides of cells
	for (int i = 1; i <= lastRealIndexX; i++)
		for (int j = 0; j <= lastRealIndexP; j++) {
			double x = getCellCenterX(i, j), p = getCellTopP(i, j);
			double wFcnVal = w_fcn(x, p);
			if (wFcnVal >= 0)
				for (int kk = 0; kk < 2; kk++)
					GG[i][j][kk] = Dx * wFcnVal * sl[i][j][kk];
			else
				for (int kk = 0; kk < 2; kk++)
					GG[i][j][kk] = Dx * wFcnVal * sl[i][j + 1][kk];
		}

	// FF fluxes: fluxes on right sides of cells
	for (int i = 0; i <= lastRealIndexX; i++)
		for (int j = 1; j <= lastRealIndexP; j++) {
			double x = getCellRightX(i, j), p = getCellCenterP(i, j);
			double uFcnVal = u_fcn(x, p);
			if (uFcnVal >= 0)
				for (int kk = 0; kk < 2; kk++)
					FF[i][j][kk] = Dp * uFcnVal * sl[i][j][kk];
			else
				for (int kk = 0; kk < 2; kk++)
					FF[i][j][kk] = Dp * uFcnVal * sl[i + 1][j][kk];
		}
}

// Calculate numerical fluxes using finite volume method
void calcFluxes_ClassFV() {
	// GG fluxes: fluxes on top sides of cells
	for (int i = 1; i <= lastRealIndexX; i++)
		for (int j = 0; j <= lastRealIndexP; j++) {
			double x = getCellCenterX(i, j), p = getCellTopP(i, j);
			for (int kk = 0; kk < 2; kk++)
				GG[i][j][kk] = Dx * w_fcn(x, p) * 0.5 * (sl[i][j][kk] + sl[i][j + 1][kk]);
		}

	// FF fluxes: fluxes on right sides of cells
	for (int i = 0; i <= lastRealIndexX; i++)
		for (int j = 1; j <= lastRealIndexP; j++) {
			double x = getCellRightX(i, j), p = getCellCenterP(i, j);
			for (int kk = 0; kk < 2; kk++)
				FF[i][j][kk] = Dp * u_fcn(x, p) * 0.5 * (sl[i][j][kk] + sl[i + 1][j][kk]);
		}
}

void (*calcFluxesPtr)();

void selectFluxFcn() {
	if (fluxMethod == 0)
		calcFluxesPtr = &calcFluxes_Godunov;
	else
		calcFluxesPtr = &calcFluxes_ClassFV;
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Time Steps
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Calculate solutions using the forward Euler method on time
void forwardEuler() {
	// Enforce the initial conditions
	enforceInitCond();

	// Select flux method
	selectFluxFcn();

	// Use the forward Euler method on time
	for (int tt = 0; tt < numTimeSteps; tt++) {
		printf("\rProgress: %1.2f%%", 100. * tt / numTimeSteps);
		double t = tt * Dt;  // Current time
		(*calcFluxesPtr)();  // Recalculate fluxes
		for (int i = 1; i < lastGhostIndexX; i++)
			for (int j = 1; j < lastGhostIndexP; j++) {
				// In each loop, forward Euler is used on the (i, j) cell
				double T = sl[i][j][0], q = sl[i][j][1], x = getCellCenterX(i, j), p = getCellCenterP(i, j);
				// RHS represents the right hand side of the formula for the forward Euler method: dT/dt = RHS
				double RHS[2] = {0, 0};
				for (int kk = 0; kk < 2; kk++)
					RHS[kk] += (GG[i][j - 1][kk] - GG[i][j][kk] + FF[i - 1][j][kk] - FF[i][j][kk]) / cellVol;
				(*addSourceFcnPtr)(RHS, T, q, x, p, t);
				for (int kk = 0; kk < 2; kk++)
					sl[i][j][kk] += Dt * RHS[kk];
			}
		(*enforceBCPtr)();  // Enforce boundary conditions
	}
	printf("\nCompleted: 100%%\n");
}

#endif /* FV_H_ */
