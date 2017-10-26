/*
 * Conditions.h
 *
 *  Created on: Oct 14, 2017
 *      Author: chuckjia
 */

#ifndef CONDITIONS_H_
#define CONDITIONS_H_
#include "Mesh.h"

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Solutions as 2D Array
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

double T_sl[numCellsX][numCellsP], q_sl[numCellsX][numCellsP];
double u_sl[numCellsX][numCellsP];
double w_sl[numCellsX][numCellsP];
double phix_sl[numCellsX][numCellsP];

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Enforce Initial Conditions
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

double (*IC_T_fcnPtr)(double x, double p, double t), (*IC_q_fcnPtr)(double x, double p, double t),
		(*IC_u_fcnPtr)(double x, double p, double t), (*IC_w_fcnPtr)(double x, double p, double t);

void enforceIC() {
	for (int i = 0; i < numCellsX; ++i)
		for (int j = 0; j < numCellsP; ++j) {
			double x = getCellCenterX(i, j), p = getCellCenterP(i, j);
			T_sl[i][j] = (*IC_T_fcnPtr)(x, p, 0);
			q_sl[i][j] = (*IC_q_fcnPtr)(x, p, 0);
			u_sl[i][j] = (*IC_u_fcnPtr)(x, p, 0);
			w_sl[i][j] = (*IC_w_fcnPtr)(x, p, 0);
		}
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Enforce Boundary Conditions
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void (*enforceBC_fcnPtr)();

// Test 1

void enforceBC_MDL1() {
	// Left boundary: Dirichlet BC
	for (int j = 1; j <= lastRealIndexP; ++j) {
		double p = getCellCenterP(1, j);
		T_sl[0][j] = 2 * exact_T_fcn_MDL1(0, p, 0) - T_sl[1][j];
		q_sl[0][j] = -q_sl[1][j];
		u_sl[0][j] = -u_sl[1][j];
	}

	// Right boundary: Neumann BC
	for (int j = 1; j <= lastRealIndexP; ++j) {
		T_sl[lastGhostIndexX][j] = T_sl[lastRealIndexX][j];
		q_sl[lastGhostIndexX][j] = q_sl[lastRealIndexX][j];
		u_sl[lastGhostIndexX][j] = u_sl[lastRealIndexX][j];
	}
}

double (*source_T_fcnPtr)(double T, double q, double u, double x, double p, double t),
		(*source_q_fcnPtr)(double T, double q, double u, double x, double p, double t),
		(*source_u_fcnPtr)(double T, double q, double u, double x, double p, double t);

void setConditions() {
	// if (modelNo == 1) {
		// Initial conditions
		IC_T_fcnPtr = &exact_T_fcn_MDL1;
		IC_q_fcnPtr = &exact_q_fcn_MDL1;
		IC_u_fcnPtr = &exact_u_fcn_MDL1;
		IC_w_fcnPtr = &exact_w_fcn_MDL1;
		// Boundary conditions
		enforceBC_fcnPtr = &enforceBC_MDL1;
		// Source functions
		source_T_fcnPtr = &source_T_fcn_MDL1;
		source_q_fcnPtr = &source_q_fcn_MDL1;
		source_u_fcnPtr = &source_u_fcn_MDL1;
}


#endif /* CONDITIONS_H_ */
