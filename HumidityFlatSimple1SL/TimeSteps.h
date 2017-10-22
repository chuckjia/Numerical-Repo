/*
 * TimeSteps.h
 *
 *  Created on: Jul 27, 2017
 *      Author: chuckjia
 */

#ifndef TIMESTEPS_H_
#define TIMESTEPS_H_
#include "Fluxes.h"

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Time Steps
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void forwardEuler() {
	enforceInitCond();
	for (int tt = 0; tt < numTimeSteps; tt++) {
		double t = tt * Dt;
		//printf("\rProgress: %1.1f%%", t / finalTime * 100);
		calcFluxes();
		for (int i = 1; i < lastGhostIndexX; i++)
			for (int j = 1; j < lastGhostIndexP; j++) {
				double T = soln[i][j], x = getCellCenterX(i, j), p = getCellCenterP(i, j);
				double RHS = source_Test1(T, x, p, t) -
						(GG_Bott[i][j] + GG_Top[i][j] + FF_Left[i][j] + FF_Right[i][j]) / cellVol;
				soln[i][j] += Dt * RHS;
			}
	}
	printf("\nCompleted: 100%%\n");
}

#endif /* TIMESTEPS_H_ */
