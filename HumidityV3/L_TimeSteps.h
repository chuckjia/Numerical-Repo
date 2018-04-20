/*
 * TimeSteps.h
 *
 *  Created on: Oct 23, 2017
 *      Author: Chuck Jia
 */

#ifndef L_TIMESTEPS_H_
#define L_TIMESTEPS_H_
#include "K_Godunov.h"

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Copies of Numerical Solutions
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

double T_sl_copy[numCellsX][numCellsP], q_sl_copy[numCellsX][numCellsP],
u_sl_copy[numCellsX][numCellsP];

double k_rk_T[numCellsX][numCellsP], k_rk_q[numCellsX][numCellsP], k_rk_u[numCellsX][numCellsP];

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Functions for RK family methods
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

double calcFluxes_OneCell(int i, int j,
		double GG[Nx + 1][Np + 1], double FF[Nx + 1][Np + 1]) {
	return GG[i][j] - GG[i][j - 1] + FF[i][j] - FF[i - 1][j];
}

void (*update_k_rk_fcnPtr)(double k_rk[numCellsX][numCellsP], int i, int j,
		double rkCoef, double kVal);

void update_k_rk_noUpdate(double k_rk[numCellsX][numCellsP], int i, int j,
		double rkCoef, double kVal) {
	// Empty
}

void update_k_rk_accum(double k_rk[numCellsX][numCellsP], int i, int j,
		double rkCoef, double kVal) {
	k_rk[i][j] += rkCoef * Dt * kVal;
}

void update_k_rk_directAssign(double k_rk[numCellsX][numCellsP], int i, int j,
		double rkCoef, double kVal) {
	k_rk[i][j] = rkCoef * Dt * kVal;
}

// Make one step using forward Euler on all real cells.
// Use data T_arr, q_arr, and u_arr
void forwardEuler_singleStep(double t, double stepSize, double T_arr[numCellsX][numCellsP],
		double q_arr[numCellsX][numCellsP], double u_arr[numCellsX][numCellsP],
		double rkCoef) {
	for (int i = 1; i <= Nx; ++i) {
		double x = getCellCenterX(i), volInv = 1. / getCellVol(i);
		for (int j = 1; j <= Np; ++j) {
			double T = T_sl[i][j], q = q_sl[i][j], u = u_sl[i][j], w = w_sl[i][j],
					p = getCellCenterP(i, j);
			double RHS;
			// Updating T
			RHS = -volInv * calcFluxes_OneCell(i, j, GG_T, FF_T) +
					(*source_T_fcnPtr)(T, q, u, w, x, p, t);
			T_sl[i][j] = T_arr[i][j] + RHS * stepSize;
			(*update_k_rk_fcnPtr)(k_rk_T, i, j, rkCoef, RHS);

			// Updating q
			RHS = -volInv * calcFluxes_OneCell(i, j, GG_q, FF_q) +
					(*source_q_fcnPtr)(T, q, u, w, x, p, t);
			q_sl[i][j] = q_arr[i][j] + RHS * stepSize;
			(*update_k_rk_fcnPtr)(k_rk_q, i, j, rkCoef, RHS);

			// Updating u
			RHS = -volInv * calcFluxes_OneCell(i, j, GG_u, FF_u) - phix_sl[i][j] +
					(*source_u_fcnPtr)(T, q, u, w, x, p, t);
			u_sl[i][j] = u_arr[i][j] + RHS * stepSize;
			(*update_k_rk_fcnPtr)(k_rk_u, i, j, rkCoef, RHS);
		}
	}
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Forward Euler Method On Time
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void forwardEuler() {
	printf("\n- Running forward Euler method on time\n");
	int prog = -1;

	// The initial condition
	enforceIC();
	update_k_rk_fcnPtr = &update_k_rk_noUpdate;

	for (int tt = 0; tt < numTimeStep; tt++) {
		// Print messages on calculation progress
		int progNew = tt * 100 / numTimeStep;
		if (progNew > prog) {
			prog = progNew; printf("\r  - Current progress: %d%%", prog); fflush(stdout);
		}

		// Numerical calculation
		double t = Dt * tt;
		(*calcFluxes)();  // Calculate numerical fluxes
		(*calc_phix_fcnPtr)();  // Calculate phi_x value at the beginning of each time step
		forwardEuler_singleStep(t, Dt, T_sl, q_sl, u_sl, 0);
		(*projU_fcnPtr)();  // Projection method on u
		(*calc_w_fcnPtr)();  // Calculate w
		(*enforceBC_fcnPtr)();  // Enforce boundary conditions

		//showL2Errors(t);
		//writeResToFileForMovie_T_test(tt + 1);
	}
	printf("\r  - Forward Euler method complete\n");
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Runge-Kutta 2
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void copySoln(double T_copy[numCellsX][numCellsP], double q_copy[numCellsX][numCellsP],
		double u_copy[numCellsX][numCellsP]) {
	for (int i = 0; i < numCellsX; ++i)
		for (int j = 0; j < numCellsP; ++j) {
			T_copy[i][j] = T_sl[i][j];
			q_copy[i][j] = q_sl[i][j];
			u_copy[i][j] = u_sl[i][j];
		}
}

void preForwardEuler() {
	enforceNeumannBC_topBD(u_sl);
	enforceNeumannBC_topBD(w_sl);
	(*calcFluxes)();  // Calculate numerical fluxes
	(*calc_phix_fcnPtr)();  // Calculate phi_x value at the beginning of each time step
}

void subExactVelocity(double t) {
	for (int i = 0; i < numCellsX; ++i) {
		double x = getCellCenterX(i);
		for (int j = 0; j < numCellsP; ++j) {
			double p = getCellCenterP(i, j);
			u_sl[i][j] = (*IC_u_fcnPtr)(x, p, t);
			w_sl[i][j] = (*IC_q_fcnPtr)(x, p, t);
		}
	}
}

void postForwardEuler() {
	(*projU_fcnPtr)();  // Projection method on u
	// (*enforceBC_fcnPtr)();  // Enforce boundary conditions
	(*calc_w_fcnPtr)();  // Calculate w
	(*enforceBC_fcnPtr)();  // Enforce boundary conditions
}

void rk2() {
	printf("\n- Running Runge-Kutta 2 method on time\n");
	int prog = -1;

	// The initial condition
	enforceIC();
	update_k_rk_fcnPtr = &update_k_rk_noUpdate;
	for (int tt = 0; tt < numTimeStep; tt++) {
		// Print messages on calculation progress
		int progNew = tt * 100 / numTimeStep;
		if (progNew > prog) {
			prog = progNew; printf("\r  - Current progress: %d%%", prog); fflush(stdout);
		}

		// Numerical calculation
		double t = Dt * tt;

		// RK2 Step 0
		copySoln(T_sl_copy, q_sl_copy, u_sl_copy);

		// RK2 Step 1
		preForwardEuler();
		forwardEuler_singleStep(t, halfDt, T_sl, q_sl, u_sl, 0);
		postForwardEuler();

		// RK2 Step 2
		preForwardEuler();
		forwardEuler_singleStep(t + halfDt, Dt, T_sl_copy, q_sl_copy, u_sl_copy, 0);
		postForwardEuler();

		(*enforceBC_fcnPtr)();

		//showL2Errors(t);
		//writeResToFileForMovie_T(tt + 1);
	}
	printf("\r  - Runge-Kutta 2 method complete\n");
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Runge-Kutta 4
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void rk4() {
	printf("\n- Running Runge-Kutta 4 method on time\n");
	int prog = -1;

	// The initial condition
	enforceIC();

	for (int tt = 0; tt < numTimeStep; tt++) {
		double t = Dt * tt;
		// Print messages on calculation progress
		int progNew = tt * 100 / numTimeStep;
		if (progNew > prog) {
			prog = progNew;
			printf("\r  - Current progress: %d%%, Time = %1.2fs", prog, t);
			fflush(stdout);
			if (calcL2err_opt)
				calcL2Norm(t);
		}

		if (!(tt % movieFrameFreq))
			printSolnToFile(tt);

		// Numerical calculation

		// RK4 Step 0
		copySoln(T_sl_copy, q_sl_copy, u_sl_copy);

		// RK4 Step 1
		update_k_rk_fcnPtr = &update_k_rk_directAssign;
		preForwardEuler();
		forwardEuler_singleStep(t, halfDt, T_sl, q_sl, u_sl, ONE_SIXTH_CONST);
		postForwardEuler();

		// RK4 Step 2
		update_k_rk_fcnPtr = &update_k_rk_accum;
		preForwardEuler();
		forwardEuler_singleStep(t + halfDt, halfDt, T_sl_copy, q_sl_copy, u_sl_copy, ONE_THIRD_CONST);
		postForwardEuler();

		// RK4 Step 3
		update_k_rk_fcnPtr = &update_k_rk_accum;
		preForwardEuler();
		forwardEuler_singleStep(t + halfDt, Dt, T_sl_copy, q_sl_copy, u_sl_copy, ONE_THIRD_CONST);
		postForwardEuler();

		// RK4 Step 4
		update_k_rk_fcnPtr = &update_k_rk_noUpdate;
		preForwardEuler();
		forwardEuler_singleStep(t + Dt, oneSixthDt, T_sl_copy, q_sl_copy, u_sl_copy, 0);
		for (int i = 1; i <= Nx; ++i)
			for (int j = 1; j <= Np; ++j) {
				T_sl[i][j] += k_rk_T[i][j];
				q_sl[i][j] += k_rk_q[i][j];
				u_sl[i][j] += k_rk_u[i][j];
			}
		postForwardEuler();

		//showL2Errors(t);
		(*aveSoln_fcnPtr)(tt + 1);

		(*enforceBC_fcnPtr)();
	}
	printSolnToFile(numTimeStep);
	printf("\r  - Runge-Kutta 4 method complete\n");
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Wrapper For Time Method
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void timeSteps() {
	clock_t start = clock();

	switch (timeMethod) {
	case 1:
		forwardEuler();
		break;
	case 2:
		rk2();
		break;
	case 4:
		rk4();
		break;
	default:
		throw "Error: Incorrect time method number!";
		return;
	}

	computationTime = (double) ((clock() - start)) / CLOCKS_PER_SEC;
	printf("\n- Calculation complete. Time used = %1.2fs.\n\n", computationTime);
}

void runTimeSteps() {
	try {
		timeSteps();
	} catch (const char* msg) {
		cerr << "\n" << msg << endl;
		exit(EXIT_FAILURE);
	}
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Set Parameters For All Time Methods
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void setTimeSteps() {
	// Empty
}


#endif /* L_TIMESTEPS_H_ */
