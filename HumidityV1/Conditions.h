/*
 * Conditions.h
 *
 *  Created on: Oct 14, 2017
 *      Author: chuckjia
 */

#ifndef CONDITIONS_H_
#define CONDITIONS_H_
#include "WPhix.h"

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
	(*projU_fcnPtr)();
	(*calc_w_fcnPtr)();
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

/**
 * LEFT
 */

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

// Cache to store the boundary values specified on the left side of boundary
double T_leftBdVal_cache[numCellsP],
q_leftBdVal_cache[numCellsP],
u_leftBdVal_cache[numCellsP];

// Enforce Dirichlet BCs on the LEFT boundary: with boundary value from cache
void enforceDirichlet_leftBD(double sl[numCellsX][numCellsP], double bdVal_cache[numCellsP]) {
	for (int j = 1; j <= Np; ++j)
		sl[0][j] = 2 * bdVal_cache[j] - sl[1][j];
}

// Enforce Neumann BCs on the left boundary
void enforceNeumann_leftBD(double sl[numCellsX][numCellsP]) {
	for (int j = 1; j <= Np; ++j)
		sl[0][j] = sl[1][j];
}

/**
 * RIGHT
 */

// Enforce Dirichlet BCs on the RIGHT boundary: with boundary value 0
void enforceDirichlet_rightBD(double sl[numCellsX][numCellsP]) {
	for (int j = 1; j <= Np; ++j)
		sl[lastGhostIndexX][j] = - sl[lastRealIndexX][j];
}

// Enforce Neumann BCs on the right boundary
void enforceNeumann_rightBD(double sl[numCellsX][numCellsP]) {
	for (int j = 1; j <= Np; ++j)
		sl[lastGhostIndexX][j] = sl[lastRealIndexX][j];
}

/**
 * BOTTOM
 */

// Enforce Dirichlet BCs on the BOTTOM boundary: with boundary value 0
void enforceDirichlet_bottBD(double sl[numCellsX][numCellsP]) {
	for (int i = 1; i <= Nx; ++i)
		sl[i][0] = - sl[i][1];
}

// Implementation from (3.39)
void enforceBC_zeroGhost_bottBD_w() {
	for (int i = 1; i <= Nx; ++i)
		w_sl[i][0] = 0;
}

/**
 * TOP
 */

void enforceNonPenetrationBC_topBD_math() {
	for (int i = 1; i <= Nx; ++i) {
		double x = getCellCenterX(i);
		w_sl[i][lastGhostIndexP] =
				pB_xDer_fcn_MDL1(x) * (u_sl[i][Np] + u_sl[i][lastGhostIndexP]) - w_sl[i][Np];
	}
}

void enforceNonPenetrationBC_topBD_numer() {
	for (int i = 1; i <= Nx; ++i) {
		double norm_x = getCellTopSideNormX(i, Np), norm_p = getCellTopSideNormP(i, Np);
		w_sl[i][lastGhostIndexP] = -norm_x * (u_sl[i][Np] + u_sl[i][lastGhostIndexP]) / norm_p
				- w_sl[i][Np];
	}
}

/* ----- ----- ----- ----- ----- -----
 * Physical Case 0: BCs
 * ----- ----- ----- ----- ----- ----- */

double helper1_fillCache_leftBdVal_MDL0(double p) {
	return sin(M_PI * p / p0_CONST);
}

double helper2_fillCache_leftBdVal_MDL0(double x) {
	return cos(TWO_PI_CONST * n_initU_coef_MDL0 / xf * x);
}

// Need to execute after enforcing initial conditions. More specifically, method only works
// after enforcing IC and applying projection on u. This is to make sure that the lambda_x
// are calculated, which will be used in the calculation.
void fillCache_leftBdVal_MDL0() {
	enforceIC();
	// Calculate lambda_x(x0)
	double pB = pB_fcn_MDL0(x0), x1 = getCellCenterX(1);
	double lambda_x_leftBdVal = lambda_x_proj[1] -
			2 * p0_CONST / (M_PI * (pB - pA))
			* (helper1_fillCache_leftBdVal_MDL0(pB) - helper1_fillCache_leftBdVal_MDL0(pA))
			* (helper2_fillCache_leftBdVal_MDL0(x1) - helper2_fillCache_leftBdVal_MDL0(x0));
	// printf("\nlambda_x on the left BD = %1.10e\n", lambda_x_leftBdVal);
	// printf("lambda_x[1] = %1.10e\n", lambda_x_proj[1]);

	for (int j = 0; j < numCellsP; ++j) {
		double p = getCellCenterP(1, j);
		// Calculate T values
		double T = init_T_fcn_MDL0(x0, p, 0);
		T_leftBdVal_cache[j] = T;
		// Calculate q values
		q_leftBdVal_cache[j] = init_q_fcn_MDL0(T, p, 0);
		// q_leftBdVal_cache[j] = qs_fcn_MDL0(T, p);
		// Calculate u values
		u_leftBdVal_cache[j] = init_u_fcn_MDL0(0, p, 0) - lambda_x_leftBdVal;
	}
}

void enforceBC_MDL0() {
	// Left boundary: Dirichlet BC
	enforceDirichlet_leftBD(T_sl, T_leftBdVal_cache);
	enforceDirichlet_leftBD(q_sl, q_leftBdVal_cache);
	enforceDirichlet_leftBD(u_sl, u_leftBdVal_cache);
	enforceNeumann_leftBD(w_sl);
	// Right boundary: Neumann BC
	enforceNeumann_rightBD(T_sl);
	enforceNeumann_rightBD(q_sl);
	enforceNeumann_rightBD(u_sl);
	enforceNeumann_rightBD(w_sl);
	// Bottom boundary: for w, Dirichlet BC with boudnary value 0
	enforceBC_zeroGhost_bottBD_w(); // Maybe not used
	// Top boundary: for u and w
	enforceNonPenetrationBC_topBD_numer();
}

/* ----- ----- ----- ----- ----- -----
 * Test Case 1 BCs
 * ----- ----- ----- ----- ----- ----- */

double leftBdVal_T_fcn_MDL1(double p) {
	return exact_T_fcn_MDL1(0, p, 0);
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
	enforceNonPenetrationBC_topBD_numer();
}

/* ----- ----- ----- ----- ----- -----
 * Test 101 BC
 * ----- ----- ----- ----- ----- ----- */

// Empty. Using Test Case 1 BC

/* ----- ----- ----- ----- ----- -----
 * Test 102 BC
 * ----- ----- ----- ----- ----- ----- */

void enforceBC_MDL102() {
	// Left boundary
	enforceDirichlet_leftBD(T_sl);
	enforceDirichlet_leftBD(q_sl);
	// Right boundary: Neumann BC
	enforceNeumann_rightBD(T_sl);
	enforceNeumann_rightBD(q_sl);
	//enforceBC_topBD_numer_MDL1();
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

// Empty for now
// Same with BCs in MDL1

/* ----- ----- ----- ----- ----- -----
 * Test 5 BC
 * ----- ----- ----- ----- ----- ----- */

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
		IC_w_fcnPtr = &zero_fcn;
		// Boundary conditions
		fillCache_leftBdVal_MDL0();
		enforceBC_fcnPtr = &enforceBC_MDL0;
		// Source functions
		source_T_fcnPtr = &source_T_fcn_MDL0;
		source_q_fcnPtr = &source_q_fcn_MDL0;
		source_u_fcnPtr = &source_u_fcn_MDL0;
		return;

	case 1:
	case 101:
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

	case 102:
		// Initial conditions
		IC_T_fcnPtr = &exact_T_fcn_MDL1;
		IC_q_fcnPtr = &exact_q_fcn_MDL1;
		IC_u_fcnPtr = &exact_u_fcn_MDL1;
		IC_w_fcnPtr = &exact_w_fcn_MDL1;
		// Boundary conditions
		enforceBC_fcnPtr = &enforceBC_MDL102;
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

	case 4:  // All are the same with MDL1, except for all T functions
		// Initial conditions
		IC_T_fcnPtr = &exact_T_fcn_MDL4;
		IC_q_fcnPtr = &exact_q_fcn_MDL1;
		IC_u_fcnPtr = &exact_u_fcn_MDL1;
		IC_w_fcnPtr = &exact_w_fcn_MDL1;
		// Boundary conditions
		enforceBC_fcnPtr = &enforceBC_MDL1;
		// Source functions
		source_T_fcnPtr = &source_T_fcn_MDL4;
		source_q_fcnPtr = &source_q_fcn_MDL1;
		source_u_fcnPtr = &source_u_fcn_MDL1;
		return;

	case 5:  // All are the same with MDL1, except for all T functions
		// Initial conditions
		IC_T_fcnPtr = &exact_T_fcn_MDL5;
		IC_q_fcnPtr = &exact_q_fcn_MDL1;
		IC_u_fcnPtr = &exact_u_fcn_MDL1;
		IC_w_fcnPtr = &exact_w_fcn_MDL1;
		// Boundary conditions
		enforceBC_fcnPtr = &enforceBC_MDL1;
		// Source functions
		source_T_fcnPtr = &source_T_fcn_MDL5;
		source_q_fcnPtr = &source_q_fcn_MDL1;
		source_u_fcnPtr = &source_u_fcn_MDL1;
		return;
	}
}

#endif /* CONDITIONS_H_ */
