/*
 * TimeSteps.h
 *
 *  Created on: Oct 23, 2017
 *      Author: chuckjia
 */

#ifndef TIMESTEPS_H_
#define TIMESTEPS_H_
#include "WPhix.h"

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Forward Euler Method On Time
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

double calcFluxOneCell(int i, int j,
		double GG[Nx + 1][Np + 1], double FF[Nx + 1][Np + 1]) {
	return GG[i][j] - GG[i][j - 1] + FF[i][j] - FF[i - 1][j];
}

void forwardEuler() {
	enforceIC();
	printf("\n- Running forward Euler method on time\n");
	for (int tt = 0; tt < numTimeSteps; tt++) {
		double t = Dt * tt;
		printf("\r  - Current progress: %1.0f%%", t / finalTime * 100);

		// Calculate phi_x value at the beginning of each time step
		calc_phix();

		// Godunov upwind
		for (int i = 1; i <= Nx; ++i) {
			double x = getCellCenterX(i), volInv = 1 / getCellVol(i);
			for (int j = 1; j <= Np; ++j) {
				double T = T_sl[i][j], q = q_sl[i][j], u = u_sl[i][j],
						p = getCellCenterP(i, j);
				double RHS;

				// Updating T
				RHS = -volInv * calcFluxOneCell(i, j, GG_T, FF_T) +
						(*source_T_fcnPtr)(T, q, u, x, p, t);
				T_sl[i][j] += RHS * Dt;
				// Updating q
				RHS = -volInv * calcFluxOneCell(i, j, GG_q, FF_q) +
						(*source_q_fcnPtr)(T, q, u, x, p, t);
				q_sl[i][j] += RHS * Dt;

				// Updating u
				RHS = -volInv * calcFluxOneCell(i, j, GG_u, FF_u) - phix_sl[i][j] +
						(*source_u_fcnPtr)(T, q, u, x, p, t);
				u_sl[i][j] += RHS * Dt;
			}
		}

		// Projection method on u
		projU();

		// Calculate w
		calc_w();

		// Enforce boundary conditions
		(*enforceBC_fcnPtr)();
	}
	printf("\r  - Forward Euler method complete\n");
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Wrapper For Time Method
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void timeSteps() {
	clock_t start = clock();

	forwardEuler();

	printf("\n- Calculation complete. Time used = %1.2fs.\n\n",
				((double) (clock() - start)) / CLOCKS_PER_SEC);
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Set Parameters For All Time Methods
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void setTimeSteps() {
	// Empty for now
}


#endif /* TIMESTEPS_H_ */
