/*
 * Conditions.h
 *
 *  Created on: Oct 14, 2017
 *      Author: chuckjia
 */

#ifndef CONDITIONS_H_
#define CONDITIONS_H_
#include "Mesh.h"

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Numerical Solutions: 2D Arrays
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

double T_sl[numCellsX][numCellsP], q_sl[numCellsX][numCellsP];
double u_sl[numCellsX][numCellsP], w_sl[numCellsX][numCellsP];
double phix_sl[numCellsX][numCellsP];

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Initial Conditions
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Function pointers: mathematical functions to calculate initial conditions
double (*IC_T_fcnPtr)(double x, double p, double t),
		(*IC_q_fcnPtr)(double x, double p, double t),
		(*IC_u_fcnPtr)(double x, double p, double t),
		(*IC_w_fcnPtr)(double x, double p, double t);  // IC for w is for testing purposes

// Enforce initial conditions at cell centers
void enforceIC() {
	// We are not setting initial conditions for w
	for (int i = 0; i < numCellsX; ++i)
		for (int j = 0; j < numCellsP; ++j) {
			double x = getCellCenterX(i, j), p = getCellCenterP(i, j);
			T_sl[i][j] = (*IC_T_fcnPtr)(x, p, 0);
			q_sl[i][j] = (*IC_q_fcnPtr)(x, p, 0);
			u_sl[i][j] = (*IC_u_fcnPtr)(x, p, 0);
		}
}

// Placeholder for w
double zeroSoln_fcn(double x, double p, double t) {
	return 0;
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Boundary Conditions
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Function pointer: function to enforce all boundary conditions in one model
void (*enforceBC_fcnPtr)();

/* ----- ----- ----- ----- ----- -----
 * Common boundary conditions
 * ----- ----- ----- ----- ----- ----- */

/*
 * In implementing all the BCs, we assume pB is sufficiently flat on the left and right domain
 * sides. This guarantees the left and right ghost cells are not distorted.
 */

// Enforce Neumann BCs on the left boundary
void enforceNeumann_leftBD(double sl[numCellsX][numCellsP]) {
	for (int j = 1; j <= Np; ++j)
		sl[0][j] = sl[1][j];
}

// Enforce Neumann BCs on the right boundary
void enforceNeumann_rightBD(double sl[numCellsX][numCellsP]) {
	for (int j = 1; j <= Np; ++j)
		sl[lastGhostIndexX][j] = sl[lastRealIndexX][j];
}

// Enforce Dirichlet BCs on the LEFT boundary: with boundary values specified by a function,
// defined by (*bdVal_fcnPtr)(p)
void enforceDirichlet_leftBD(double sl[numCellsX][numCellsP], double (&bdVal_fcn)(double)) {
	for (int j = 1; j <= Np; ++j) {
		// This interpolation on p only works for pB() that are sufficiently flat on the left
		// side of domain
		double p = 0.5 * (getCellBottLeftP(1, j) + getCellTopLeftP(1, j));
		sl[0][j] = 2 * bdVal_fcn(p) - sl[1][j];
	}
}

// Enforce Dirichlet BCs on the LEFT boundary: with constant boundary value, specified by bdVal
void enforceDirichlet_leftBD(double sl[numCellsX][numCellsP], double bdVal) {
	for (int j = 1; j <= Np; ++j)
		sl[0][j] = 2 * bdVal - sl[1][j];
}

// Enforce Dirichlet BCs on the LEFT boundary: with boundary value 0
void enforceDirichlet_leftBD(double sl[numCellsX][numCellsP]) {
	for (int j = 1; j <= Np; ++j)
		sl[0][j] = - sl[1][j];
}

// Enforce Dirichlet BCs on the RIGHT boundary: with boundary value 0
void enforceDirichlet_rightBD(double sl[numCellsX][numCellsP]) {
	for (int j = 1; j <= Np; ++j)
		sl[lastGhostIndexX][j] = - sl[lastRealIndexX][j];
}

// Enforce Dirichlet BCs on the BOTTOM boundary: with boundary value 0
void enforceDirichlet_bottBD(double sl[numCellsX][numCellsP]) {
	for (int i = 1; i <= Nx; ++i)
		sl[i][0] = - sl[i][1];
}

/* ----- ----- ----- ----- ----- -----
 * Physical Case 0 BCs
 * ----- ----- ----- ----- ----- ----- */

double leftBdVal_T_fcn_MDL0(double p) {
	return init_T_fcn_MDL0(0, p, 0);
}

double leftBdVal_q_fcn_MDL0(double p) {
	return qs_fcn_MDL0(init_T_fcn_MDL0(0, p, 0), p);
}

// Not complete
double leftBdVal_u_fcn_MDL0(double p) {
	return 0;
}

void enforceBC_MDL0() {
	// Left boundary: Dirichlet BC
	enforceDirichlet_leftBD(T_sl, leftBdVal_T_fcn_MDL0);
	enforceDirichlet_leftBD(q_sl, leftBdVal_q_fcn_MDL0);
	enforceDirichlet_leftBD(u_sl, leftBdVal_u_fcn_MDL0);
	enforceDirichlet_leftBD(w_sl);  // Maybe not used: WRONG!
	// Right boundary: Neumann BC
	enforceNeumann_rightBD(T_sl);
	enforceNeumann_rightBD(q_sl);
	enforceNeumann_rightBD(u_sl);
	enforceNeumann_rightBD(w_sl);
	// Bottom boundary: for w, Dirichlet BC with boudnary value 0
	enforceDirichlet_bottBD(w_sl);  // Maybe not used
	// Top boundary: for u and w
}

/* ----- ----- ----- ----- ----- -----
 * Test Case 1 BCs
 * ----- ----- ----- ----- ----- ----- */

double leftBdVal_T_fcn_MDL1(double p) {
	return exact_T_fcn_MDL1(0, p, 0);
}

void enforceBC_topBD_math_MDL1() {
	for (int i = 1; i <= Nx; ++i) {
		double x = getCellCenterX(i);
		w_sl[i][lastGhostIndexP] =
				pB_xDer_fcn_MDL1(x) * (u_sl[i][Np] + u_sl[i][lastGhostIndexP]) - w_sl[i][Np];
	}
}

void enforceBC_topBD_numer_MDL1() {
	for (int i = 1; i <= Nx; ++i) {
		double norm_x = getCellTopSideNormX(i, Np), norm_p = getCellTopSideNormP(i, Np);
		w_sl[i][lastGhostIndexP] = -norm_x * (u_sl[i][Np] + u_sl[i][lastGhostIndexP])
														 / norm_p - w_sl[i][Np];
	}
}

// Implementation from (3.39)
void enforceBC_bottBD_w_MDL1() {
	for (int i = 1; i <= Nx; ++i)
		w_sl[i][0] = 0;
}

void enforceBC_MDL1() {
	// Left boundary
	// enforceDirichlet_leftBD(T_sl, leftBdVal_T_fcn_MDL1);
	enforceDirichlet_leftBD(T_sl);
	enforceDirichlet_leftBD(q_sl);
	enforceDirichlet_leftBD(u_sl);
	enforceNeumann_leftBD(w_sl);
	// Right boundary: Neumann BC
	enforceNeumann_rightBD(T_sl);
	enforceNeumann_rightBD(q_sl);
	enforceNeumann_rightBD(u_sl);
	enforceNeumann_rightBD(w_sl);
	// Bottom boundary: for w, Dirichlet BC with boudnary value 0
	// enforceBC_bottBD_w_MDL1();  // Maybe not used. From (3.39)
	// Top boundary: for u and w
	//enforceBC_topBD_math_MDL1();
	enforceBC_topBD_numer_MDL1();
}

/* ----- ----- ----- ----- ----- -----
 * Test 2 BC
 * ----- ----- ----- ----- ----- ----- */

void enforceBC_MDL2() {
	// Left boundary: Dirichlet BC
	enforceDirichlet_leftBD(T_sl, 1);
	enforceDirichlet_leftBD(q_sl, 1);
	enforceDirichlet_leftBD(u_sl, 1);
	// Right boundary: Neumann BC
	enforceNeumann_rightBD(T_sl);
	enforceNeumann_rightBD(q_sl);
	enforceNeumann_rightBD(u_sl);
}

/* ----- ----- ----- ----- ----- -----
 * Test 3 BC
 * ----- ----- ----- ----- ----- ----- */

void enforceBC_MDL3() {
	// Left boundary: Dirichlet BC
	enforceDirichlet_leftBD(T_sl);
	enforceDirichlet_leftBD(q_sl);
	// Right boundary: Dirichlet BC
	enforceDirichlet_rightBD(T_sl);
	enforceDirichlet_rightBD(q_sl);
}

/* ----- ----- ----- ----- ----- -----
 * Test 4 BC
 * ----- ----- ----- ----- ----- ----- */

void enforceBC_MDL4() {
	// Left boundary
	enforceDirichlet_leftBD(T_sl);
	enforceDirichlet_leftBD(q_sl);
	enforceDirichlet_leftBD(u_sl);
	enforceNeumann_leftBD(w_sl);
	// Right boundary: Neumann BC
	enforceNeumann_rightBD(T_sl);
	enforceNeumann_rightBD(q_sl);
	enforceNeumann_rightBD(u_sl);
	enforceNeumann_rightBD(w_sl);
	// Bottom boundary: for w, Dirichlet BC with boudnary value 0
	// enforceBC_bottBD_w_MDL1();  // Maybe not used. From (3.39)
	// Top boundary: for u and w
	//enforceBC_topBD_math_MDL1();
	enforceBC_topBD_numer_MDL1();
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Source Functions
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Function pointers: mathematical functions to calculate source term values
double (*source_T_fcnPtr)(double T, double q, double u, double w, double x, double p, double t),
		(*source_q_fcnPtr)(double T, double q, double u, double w, double x, double p, double t),
		(*source_u_fcnPtr)(double T, double q, double u, double w, double x, double p, double t);

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Set All Conditions: ICs, BCs and Source Functions
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void setConditions() {
	switch (modelNo) {
	case 0:
		// Initial conditions
		IC_T_fcnPtr = &init_T_fcn_MDL0;
		IC_q_fcnPtr = &init_q_fcn_MDL0;
		IC_u_fcnPtr = &init_u_fcn_MDL0;
		IC_w_fcnPtr = &zeroSoln_fcn;
		// Boundary conditions
		enforceBC_fcnPtr = &enforceBC_MDL0;
		// Source functions
		source_T_fcnPtr = &source_T_fcn_MDL0;
		source_q_fcnPtr = &source_q_fcn_MDL0;
		source_u_fcnPtr = &source_u_fcn_MDL0;
		return;

	case 1:
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
		return;

	case 2:
		// Initial conditions
		IC_T_fcnPtr = &exact_T_fcn_MDL2;
		IC_q_fcnPtr = &exact_q_fcn_MDL2;
		IC_u_fcnPtr = &exact_u_fcn_MDL2;
		IC_w_fcnPtr = &exact_w_fcn_MDL2;
		// Boundary conditions
		enforceBC_fcnPtr = &enforceBC_MDL2;
		// Source functions
		source_T_fcnPtr = &source_T_fcn_MDL2;
		source_q_fcnPtr = &source_q_fcn_MDL2;
		source_u_fcnPtr = &source_u_fcn_MDL2;
		return;

	case 3:
		// Initial conditions
		IC_T_fcnPtr = &exact_T_fcn_MDL3;
		IC_q_fcnPtr = &exact_q_fcn_MDL3;
		IC_u_fcnPtr = &exact_u_fcn_MDL3;
		IC_w_fcnPtr = &exact_w_fcn_MDL3;
		// Boundary conditions
		enforceBC_fcnPtr = &enforceBC_MDL3;
		// Source functions
		source_T_fcnPtr = &source_T_fcn_MDL3;
		source_q_fcnPtr = &source_q_fcn_MDL3;
		source_u_fcnPtr = &source_u_fcn_MDL3;
		return;

	case 4:
		// Initial conditions
		IC_T_fcnPtr = &exact_T_fcn_MDL4;
		IC_q_fcnPtr = &exact_q_fcn_MDL4;
		IC_u_fcnPtr = &exact_u_fcn_MDL4;
		IC_w_fcnPtr = &exact_w_fcn_MDL4;
		// Boundary conditions
		enforceBC_fcnPtr = &enforceBC_MDL4;
		// Source functions
		source_T_fcnPtr = &source_T_fcn_MDL4;
		source_q_fcnPtr = &source_q_fcn_MDL4;
		source_u_fcnPtr = &source_u_fcn_MDL4;
		return;
	}
}

#endif /* CONDITIONS_H_ */
