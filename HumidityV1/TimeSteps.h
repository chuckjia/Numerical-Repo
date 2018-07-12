/*
 * TimeSteps.h
 *
 *  Created on: Oct 23, 2017
 *      Author: Chuck Jia
 */

#ifndef TIMESTEPS_H_
#define TIMESTEPS_H_
#include "Godunov.h"

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Copies of Numerical Solutions
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

double T_copy_[numCellX][numCellP], q_copy_[numCellX][numCellP],
u_copy_[numCellX][numCellP];

double k_rk_T_[numCellX][numCellP], k_rk_q_[numCellX][numCellP], k_rk_u_[numCellX][numCellP];

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Functions for RK family methods
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

double calcFluxes_OneCell(int i, int j,
		double GG[Nx + 1][Np + 1], double FF[Nx + 1][Np + 1]) {
	return GG[i][j] - GG[i][j - 1] + FF[i][j] - FF[i - 1][j];
}

void (*update_k_RK_fptr)(double k_rk[numCellX][numCellP], int i, int j,
		double rkCoef, double kVal);

void update_k_RK_noUpdate(double k_rk[numCellX][numCellP], int i, int j,
		double rkCoef, double kVal) {
	// Empty
}

void update_k_RK_accum(double k_rk[numCellX][numCellP], int i, int j,
		double rkCoef, double kVal) {
	k_rk[i][j] += rkCoef * Dt * kVal;
}

void update_k_RK_directAssign(double k_rk[numCellX][numCellP], int i, int j,
		double rkCoef, double kVal) {
	k_rk[i][j] = rkCoef * Dt * kVal;
}

// Make one step using forward Euler on all real cells.
// Use data T_arr, q_arr, and u_arr
void forwardEuler_singleStep(double t, double stepSize, double T_arr[numCellX][numCellP],
		double q_arr[numCellX][numCellP], double u_arr[numCellX][numCellP],
		double rkCoef) {
	for (int i = 1; i <= Nx; ++i) {
		double x = getCellCenterX(i), volInv = 1. / getCellVol(i);
		for (int j = 1; j <= Np; ++j) {
			double T = T_[i][j], q = q_[i][j], u = u_[i][j], w = w_[i][j],
					p = getCellCenterP(i, j);
			double RHS;
			// Updating T
			RHS = -volInv * calcFluxes_OneCell(i, j, GG_T_, FF_T_) +
					(*source_T_fcnPtr)(T, q, u, w, x, p, t);
			T_[i][j] = T_arr[i][j] + RHS * stepSize;
			(*update_k_RK_fptr)(k_rk_T_, i, j, rkCoef, RHS);

			// Updating q
			RHS = -volInv * calcFluxes_OneCell(i, j, GG_q_, FF_q_) +
					(*source_q_fcnPtr)(T, q, u, w, x, p, t);
			q_[i][j] = q_arr[i][j] + RHS * stepSize;
			(*update_k_RK_fptr)(k_rk_q_, i, j, rkCoef, RHS);

			// Updating u
			RHS = -volInv * calcFluxes_OneCell(i, j, GG_u_, FF_u_) - phix_[i][j] +
					(*source_u_fcnPtr)(T, q, u, w, x, p, t);
			u_[i][j] = u_arr[i][j] + RHS * stepSize;
			(*update_k_RK_fptr)(k_rk_u_, i, j, rkCoef, RHS);
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
	update_k_RK_fptr = &update_k_RK_noUpdate;

	for (int tt = 0; tt < numTimeStep; tt++) {
		// Print messages on calculation progress
		int progNew = tt * 100 / numTimeStep;
		if (progNew > prog) {
			prog = progNew; printf("\r  - Current progress: %d%%", prog); fflush(stdout);
		}

		// Numerical calculation
		double t = Dt * tt;
		(*calcFluxes)();  // Calculate numerical fluxes
		(*calcPhix_fptr)();  // Calculate phi_x value at the beginning of each time step
		forwardEuler_singleStep(t, Dt, T_, q_, u_, 0);
		(*projU_fptr)();  // Projection method on u
		(*calcW_fptr)();  // Calculate w
		(*enforceBC_fptr)();  // Enforce boundary conditions

		//showL2Errors(t);
		//writeResToFileForMovie_T_test(tt + 1);
	}
	printf("\r  - Forward Euler method complete\n");
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Runge-Kutta 2
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void copySoln(double T_copy[numCellX][numCellP], double q_copy[numCellX][numCellP],
		double u_copy[numCellX][numCellP]) {
	for (int i = 0; i < numCellX; ++i)
		for (int j = 0; j < numCellP; ++j) {
			T_copy[i][j] = T_[i][j];
			q_copy[i][j] = q_[i][j];
			u_copy[i][j] = u_[i][j];
		}
}

void preForwardEuler() {
	(*calcFluxes)();  // Calculate numerical fluxes
	(*calcPhix_fptr)();  // Calculate phi_x value at the beginning of each time step
}

void subExactVelocity(double t) {
	for (int i = 0; i < numCellX; ++i) {
		double x = getCellCenterX(i);
		for (int j = 0; j < numCellP; ++j) {
			double p = getCellCenterP(i, j);
			u_[i][j] = (*initU_fptr)(x, p, t);
			w_[i][j] = (*initQ_fptr)(x, p, t);
		}
	}
}

void postForwardEuler() {
	(*projU_fptr)();  // Projection method on u
	(*enforceBC_fptr)();  // Enforce boundary conditions
	(*calcW_fptr)();  // Calculate w
	(*enforceBC_fptr)();  // Enforce boundary conditions
}

void rk2() {
	printf("\n- Running Runge-Kutta 2 method on time\n");
	int prog = -1;

	// The initial condition
	enforceIC();
	update_k_RK_fptr = &update_k_RK_noUpdate;
	for (int tt = 0; tt < numTimeStep; tt++) {
		// Print messages on calculation progress
		int progNew = tt * 100 / numTimeStep;
		if (progNew > prog) {
			prog = progNew; printf("\r  - Current progress: %d%%", prog); fflush(stdout);
		}

		// Numerical calculation
		double t = Dt * tt;

		// RK2 Step 0
		copySoln(T_copy_, q_copy_, u_copy_);

		// RK2 Step 1
		preForwardEuler();
		forwardEuler_singleStep(t, halfDt, T_, q_, u_, 0);
		postForwardEuler();

		// RK2 Step 2
		preForwardEuler();
		forwardEuler_singleStep(t + halfDt, Dt, T_copy_, q_copy_, u_copy_, 0);
		postForwardEuler();

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
			if (_calcL2err_)
				calcL2Norm(t);
		}

		writeMovie_soln(tt);

		// Numerical calculation

		// RK4 Step 0
		copySoln(T_copy_, q_copy_, u_copy_);

		// RK4 Step 1
		update_k_RK_fptr = &update_k_RK_directAssign;
		preForwardEuler();
		forwardEuler_singleStep(t, halfDt, T_, q_, u_, ONE_SIXTH);
		postForwardEuler();

		// RK4 Step 2
		update_k_RK_fptr = &update_k_RK_accum;
		preForwardEuler();
		forwardEuler_singleStep(t + halfDt, halfDt, T_copy_, q_copy_, u_copy_, ONE_THIRD);
		postForwardEuler();

		// RK4 Step 3
		update_k_RK_fptr = &update_k_RK_accum;
		preForwardEuler();
		forwardEuler_singleStep(t + halfDt, Dt, T_copy_, q_copy_, u_copy_, ONE_THIRD);
		postForwardEuler();

		// RK4 Step 4
		update_k_RK_fptr = &update_k_RK_noUpdate;
		preForwardEuler();
		forwardEuler_singleStep(t + Dt, oneSixthDt, T_copy_, q_copy_, u_copy_, 0);
		for (int i = 1; i <= Nx; ++i)
			for (int j = 1; j <= Np; ++j) {
				T_[i][j] += k_rk_T_[i][j];
				q_[i][j] += k_rk_q_[i][j];
				u_[i][j] += k_rk_u_[i][j];
			}
		postForwardEuler();

		//showL2Errors(t);
		//writeResToFileForMovie_T_test(tt + 1);
		(*aveSoln_fptr)(tt + 1);
	}
	writeMovie_soln(numTimeStep);
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


#endif /* TIMESTEPS_H_ */
