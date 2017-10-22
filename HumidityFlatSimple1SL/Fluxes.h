/*
 * Fluxes.h
 *
 *  Created on: Jul 18, 2017
 *      Author: chuckjia
 */

#ifndef FLUXES_H_
#define FLUXES_H_
#include "Models.h"

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Fluxes
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

double FF_Left[numCellsX][numCellsP], FF_Right[numCellsX][numCellsP];
double GG_Bott[numCellsX][numCellsP], GG_Top[numCellsX][numCellsP];

void calcFluxes_Dirichlet() {
	for (int i = 1; i < lastGhostIndexX; i++)
		for (int j = 1; j < lastGhostIndexP; j++) {
			double x, p;
			// GG fluxes: bottom
			x = getCellCenterX(i, j);
			p = getCellBottP(i, j);
			if (j == 1)
				GG_Bott[i][j] = -Dx * w_fcn(x, p) * boundaryVal;
			else
				GG_Bott[i][j] = -Dx * w_fcn(x, p) * 0.5 * (soln[i][j - 1] + soln[i][j]);

			// GG fluxes: top
			p = getCellTopP(i, j);
			if (j == lastRealIndexP)
				GG_Top[i][j] = Dx * w_fcn(x, p) * boundaryVal;
			else
				GG_Top[i][j] = Dx * w_fcn(x, p) * 0.5 * (soln[i][j] + soln[i][j + 1]);

			// FF fluxes: left
			p = getCellCenterP(i, j);
			x = getCellLeftX(i, j);
			if (i == 1)
				FF_Left[i][j] = -Dp * u_fcn(x, p) * boundaryVal;
			else
				FF_Left[i][j] = -Dp * u_fcn(x, p) * 0.5 * (soln[i - 1][j] + soln[i][j]);

			// FF fluxes: right
			x = getCellRightX(i, j);
			if (i == lastRealIndexX)
				FF_Right[i][j] = Dp * u_fcn(x, p) * boundaryVal;
			else
				FF_Right[i][j] = Dp * u_fcn(x, p) * 0.5 * (soln[i][j] + soln[i + 1][j]);
		}
}

void calcFluxes() {
	return calcFluxes_Dirichlet();
}

#endif /* FLUXES_H_ */
