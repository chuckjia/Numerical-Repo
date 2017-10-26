/*
 * TimeSteps.h
 *
 *  Created on: Oct 23, 2017
 *      Author: chuckjia
 */

#ifndef TIMESTEPS_H_
#define TIMESTEPS_H_
#include "WPhix.h"

double calcFluxOneCell(int i, int j,
		double GG[Nx + 1][Np + 1], double FF[Nx + 1][Np + 1]) {
	return GG[i][j] - GG[i][j - 1] + FF[i][j] - FF[i][j - 1];
}

void forwardEuler() {
	for (int tt = 0; tt < numTimeSteps; tt++) {
		double t = Dt * tt;
		for (int i = 1; i <= lastRealIndexX; ++i) {
			double x = getCellCenterX(i);
			for (int j = 1; j <= lastRealIndexP; ++j) {
				double T = T_sl[i][j], q = q_sl[i][j], u = u_sl[i][j],
						p = getCellCenterP(i, j), volInv = 1 / getCellVol(i, j);
				double R;

				// Updating T
				R = -volInv * calcFluxOneCell(i, j, GG_T, FF_T)
				+ (*source_T_fcnPtr)(T, q, u, x, p, t);
				T_sl[i][j] += R * Dt;
				// Updating q
				R = -volInv * calcFluxOneCell(i, j, GG_q, FF_q)
				+ (*source_q_fcnPtr)(T, q, u, x, p, t);
				q_sl[i][j] += R * Dt;

				// Updating u
				R = -volInv * calcFluxOneCell(i, j, GG_u, FF_u)
				+ (*source_u_fcnPtr)(T, q, u, x, p, t) - phix_sl[i][j];
				u_sl[i][j] += R * Dt;
				// Projection method
				projU();
			}
		}

	}
}

void timeSteps() {
	forwardEuler();
}

void setTimeSteps() {
	// Empty for now
}


#endif /* TIMESTEPS_H_ */
