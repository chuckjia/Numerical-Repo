/*
 * Fluxes.h
 *
 *  Created on: Jul 14, 2017
 *      Author: chuckjia
 */

#ifndef FLUXES_H_
#define FLUXES_H_
#include "Conditions.h"

/*
 * Initialize the fluxes
 */
double FF[numCellsX][numCellsP][2];
double GG[numCellsX][numCellsP][2];

double dirichletBoundaryVal;

const int secLastIndexX = Nx;
const int secLastIndexP = Np;

void calcFluxes_ClassFV_Dirichlet() {
	for (int i = 1; i < secLastIndexX; i++)
		for (int j = 1; j < secLastIndexP; j++)
			for (int ii = 0; ii < 2; ii++) {
				double x = getCellRightX(i, j), p = getCellCenterP(i, j);
				FF[i][j][ii] = 0.5 * (soln[i][j][ii] + soln[i + 1][j][ii]) * u_fcn(x, p) * Dp;
				x = getCellCenterX(i, j), p = getCellTopP(i, j);
				GG[i][j][ii] = 0.5 * (soln[i][j][ii] + soln[i][j + 1][ii]) * omega_fcn(x, p) * Dx;
			}
	for (int j = 0; j < lastGhostP; j++)
		for (int ii = 0; ii < 2; ii++) {
			FF[0][j][ii] = dirichletBoundaryVal;
			FF[secLastIndexX][j][ii] = dirichletBoundaryVal;
		}
	for (int i = 0; i < lastGhostX; i++)
		for (int ii = 0; ii < 2; ii++) {
			GG[i][0][ii] = dirichletBoundaryVal;
			GG[i][secLastIndexP][ii] = dirichletBoundaryVal;
		}
}

/*
 * Calculate and add to the answer the RHS of the RK method
 */
void addRHSofRK_FV(double ans[2], int i, int j) {
	double cellVolVal = getCellVol(i, j);
	for (int ii = 0; ii < 2; ii++)
		ans[ii] += - (GG[i][j][ii] - GG[i][j - 1][ii] +
				FF[i][j][ii] - FF[i - 1][j][ii]) / cellVolVal;
}

#endif /* FLUXES_H_ */
