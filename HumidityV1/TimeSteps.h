/*
 * TimeSteps.h
 *
 *  Created on: Oct 23, 2017
 *      Author: Chuck Jia
 */

#ifndef TIMESTEPS_H_
#define TIMESTEPS_H_
#include "WPhix.h"

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Forward Euler Method On Time
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

double calcFluxes_OneCell(int i, int j,
		double GG[Nx + 1][Np + 1], double FF[Nx + 1][Np + 1]) {
	return GG[i][j] - GG[i][j - 1] + FF[i][j] - FF[i - 1][j];
}

// Make one step using forward Euler on all real cells.
// Use data T_arr, q_arr, and u_arr
void forwardEuler_singleStep(double t, double T_arr[numCellsX][numCellsP],
		double q_arr[numCellsX][numCellsP], double u_arr[numCellsX][numCellsP]) {
	for (int i = 1; i <= Nx; ++i) {
		double x = getCellCenterX(i), volInv = 1. / getCellVol(i);
		for (int j = 1; j <= Np; ++j) {
			double T = T_sl[i][j], q = q_sl[i][j], u = u_sl[i][j],
					p = getCellCenterP(i, j);
			double RHS;
			// Updating T
			RHS = -volInv * calcFluxes_OneCell(i, j, GG_T, FF_T) +
					(*source_T_fcnPtr)(T, q, u, x, p, t);
			T_sl[i][j] = T_arr[i][j] + RHS * Dt;
			// Updating q
			RHS = -volInv * calcFluxes_OneCell(i, j, GG_q, FF_q) +
					(*source_q_fcnPtr)(T, q, u, x, p, t);
			q_sl[i][j] = q_arr[i][j] + RHS * Dt;
			// Updating u
			RHS = -volInv * calcFluxes_OneCell(i, j, GG_u, FF_u) - phix_sl[i][j] +
					(*source_u_fcnPtr)(T, q, u, x, p, t);
			u_sl[i][j] = u_arr[i][j] + RHS * Dt;
		}
	}
}

void forwardEuler() {
	printf("\n- Running forward Euler method on time\n");
	int prog = -1;

	// The initial condition
	enforceIC();
	for (int tt = 0; tt < numTimeSteps; tt++) {
		// Print messages on calculation progress
		int progNew = tt * 100 / numTimeSteps;
		if (progNew > prog) {
			prog = progNew; printf("\r  - Current progress: %d%%", prog); fflush(stdout);
		}

		// Numerical calculation
		double t = Dt * tt;
		(*calcFluxes)();  // Calculate numerical fluxes
		(*calc_phix_fcnPtr)();  // Calculate phi_x value at the beginning of each time step
		forwardEuler_singleStep(t, T_sl, q_sl, u_sl);
		(*projU_fcnPtr)();  // Projection method on u
		(*calc_w_fcnPtr)();  // Calculate w
		(*enforceBC_fcnPtr)();  // Enforce boundary conditions

		//showL2Errors(t);
		//writeResToFileForMovie_T(tt + 1);
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
