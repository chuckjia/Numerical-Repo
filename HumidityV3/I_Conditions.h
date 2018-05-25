/*
 * Conditions.h
 *
 *  Created on: Oct 14, 2017
 *      Author: chuckjia
 */

#ifndef I_CONDITIONS_H_
#define I_CONDITIONS_H_
#include "H_WPhix.h"

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Initial Conditions
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Function pointers: mathematical functions to calculate initial conditions
double (*IC_T_fcnPtr)(double x, double p, double t),
(*IC_q_fcnPtr)(double x, double p, double t),
(*IC_u_fcnPtr)(double x, double p, double t),
(*IC_w_fcnPtr)(double x, double p, double t);  // IC for w is for testing purposes

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
void enforceDirichlet_leftBD(double sl[numCellX][numCellP], double (&bdVal_fcn)(double)) {
	for (int j = 1; j <= Np; ++j) {
		// This interpolation on p only works for pB() that are sufficiently flat on the left
		// side of domain
		double p = 0.5 * (getCellBottLeftP(1, j) + getCellTopLeftP(1, j));
		sl[0][j] = 2 * bdVal_fcn(p) - sl[1][j];
	}
}

// Enforce Dirichlet BCs on the LEFT boundary: with constant boundary value, specified by bdVal
void enforceDirichlet_leftBD(double sl[numCellX][numCellP], double bdVal) {
	for (int j = 1; j <= Np; ++j)
		sl[0][j] = 2 * bdVal - sl[1][j];
}

// Enforce Dirichlet BCs on the LEFT boundary: with boundary value 0
void enforceDirichlet_leftBD(double sl[numCellX][numCellP]) {
	for (int j = 1; j <= Np; ++j)
		sl[0][j] = - sl[1][j];
}

// Cache to store the boundary values specified on the left side of boundary
double T_leftBdVal_cache[numCellP],
q_leftBdVal_cache[numCellP],
u_leftBdVal_cache[numCellP];

// Enforce Dirichlet BCs on the LEFT boundary: with boundary value from cache
void enforceDirichlet_leftBD(double sl[numCellX][numCellP], double bdVal_cache[numCellP]) {
	for (int j = 1; j <= Np; ++j)
		sl[0][j] = 2 * bdVal_cache[j] - sl[1][j];
}

// Enforce Neumann BCs on the left boundary
void enforceNeumann_leftBD(double sl[numCellX][numCellP]) {
	for (int j = 1; j <= Np; ++j)
		sl[0][j] = sl[1][j];
}

/**
 * RIGHT
 */

// Enforce Dirichlet BCs on the RIGHT boundary: with boundary value 0
void enforceDirichlet_rightBD(double sl[numCellX][numCellP]) {
	for (int j = 1; j <= Np; ++j)
		sl[lastGhostIndexX][j] = - sl[lastRealIndexX][j];
}

// Enforce Neumann BCs on the right boundary
void enforceNeumann_rightBD(double sl[numCellX][numCellP]) {
	for (int j = 1; j <= Np; ++j)
		sl[lastGhostIndexX][j] = sl[lastRealIndexX][j];
}

/**
 * BOTTOM
 */

// Enforce Dirichlet BCs on the BOTTOM boundary: with boundary value 0
void enforceDirichlet_bottBD(double sl[numCellX][numCellP]) {
	for (int i = 1; i <= Nx; ++i)
		sl[i][0] = - sl[i][1];
}

// Implementation from (3.39)
void enforceBC_zeroGhost_bottBD_w() {
	for (int i = 1; i <= Nx; ++i)
		w_[i][0] = 0;
}

/**
 * TOP
 */

void enforceNonPenetrationBC_topBD_math() {
	for (int i = 1; i <= Nx; ++i) {
		double x = getCellCenterX(i);
		w_[i][lastGhostIndexP] =
				pBxDer_fcn_MDL1(x) * (u_[i][Np] + u_[i][lastGhostIndexP]) - w_[i][Np];
	}
}

void enforceNonPenetrationBC_topBD_numer() {
	for (int i = 1; i <= Nx; ++i) {
		double norm_x = getCellTopSideNormVecX(i, Np), norm_p = getCellTopSideNormVecP(i, Np);
		w_[i][lastGhostIndexP] = -norm_x * (u_[i][Np] + u_[i][lastGhostIndexP]) / norm_p
				- w_[i][Np];
	}
}

void enforceNeumannBC_topBD(double sol[numCellX][numCellP]) {
	for (int i = 1; i <= Nx; ++i)
		sol[i][lastGhostIndexP] = sol[i][Np];
}

// Enforce initial conditions at cell centers
void enforceIC() {
	// We are not setting initial conditions for w
	for (int i = 0; i < numCellX; ++i)
		for (int j = 0; j < numCellP; ++j) {
			double x = getCellCenterX(i, j), p = getCellCenterP(i, j);
			T_[i][j] = (*IC_T_fcnPtr)(x, p, 0);
			q_[i][j] = (*IC_q_fcnPtr)(x, p, 0);
			u_[i][j] = (*IC_u_fcnPtr)(x, p, 0);
		}
	(*projU_fcnPtr)();
	(*calc_w_fcnPtr)();
	enforceBC_zeroGhost_bottBD_w(); // Maybe not used
	// enforceNonPenetrationBC_topBD_numer();
	enforceNeumannBC_topBD(u_);
	enforceNeumannBC_topBD(w_);
}

/* ----- ----- ----- ----- ----- -----
 * Physical Case 0: BCs
 * ----- ----- ----- ----- ----- ----- */

double helper1_fillCache_leftBdVal_MDL0(double p) {
	return sin(M_PI * p / p0_CONST);
}

double helper2_fillCache_leftBdVal_MDL0(double x) {
	return cos(TWO_PI * _n_initU_MDL0 / xf * x);
}

// Need to execute after enforcing initial conditions. More specifically, method only works
// after enforcing IC and applying projection on u. This is to make sure that the lambda_x
// are calculated, which will be used in the calculation.
void fillCache_leftBdVal_MDL0() {
	// The following is to project u for the left boundary values
	// enforceIC();
	// Calculate lambda_x(x0)
	// double pB = pB_fcn_MDL0(x0), x1 = getCellCenterX(1);
	double lambda_x_leftBdVal = 0; // lambda_x_proj[1] -
	//2 * p0_CONST / (M_PI * (pB - pA))
	//* (helper1_fillCache_leftBdVal_MDL0(pB) - helper1_fillCache_leftBdVal_MDL0(pA))
	//* (helper2_fillCache_leftBdVal_MDL0(x1) - helper2_fillCache_leftBdVal_MDL0(x0));
	//printf("\nlambda_x on the left BD = %1.10e\n", lambda_x_leftBdVal);
	//printf("lambda_x[1] = %1.10e\n", lambda_x_proj[1]);

	for (int j = 0; j < numCellP; ++j) {
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
	enforceDirichlet_leftBD(T_, T_leftBdVal_cache);
	enforceDirichlet_leftBD(q_, q_leftBdVal_cache);
	enforceDirichlet_leftBD(u_, u_leftBdVal_cache);
	enforceNeumann_leftBD(w_);
	// Right boundary: Neumann BC
	enforceNeumann_rightBD(T_);
	enforceNeumann_rightBD(q_);
	enforceNeumann_rightBD(u_);
	enforceNeumann_rightBD(w_);
	// Bottom boundary: for w, Dirichlet BC with boudnary value 0
	enforceBC_zeroGhost_bottBD_w(); // Maybe not used
	// Top boundary: for u and w
	// enforceNonPenetrationBC_topBD_numer();
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
	enforceDirichlet_leftBD(T_);
	enforceDirichlet_leftBD(q_);
	enforceDirichlet_leftBD(u_);
	enforceNeumann_leftBD(w_);
	// Right boundary: Neumann BC
	enforceNeumann_rightBD(T_);
	enforceNeumann_rightBD(q_);
	enforceNeumann_rightBD(u_);
	enforceNeumann_rightBD(w_);
	// Bottom boundary: for w, Dirichlet BC with boudnary value 0
	// enforceBC_bottBD_w_MDL1();  // Maybe not used. From (3.39)
	// Top boundary: for u and w
	//enforceBC_topBD_math_MDL1();
	enforceNonPenetrationBC_topBD_numer();
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
		// Initial conditions
		IC_T_fcnPtr = &exact_T_fcn_MDL1;
		IC_q_fcnPtr = &exact_q_fcn_MDL1;
		IC_u_fcnPtr = &exact_U_fcn_MDL1;
		IC_w_fcnPtr = &exact_W_fcn_MDL1;
		// Boundary conditions
		enforceBC_fcnPtr = &enforceBC_MDL1;
		// Source functions
		source_T_fcnPtr = &source_T_fcn_MDL1;
		source_q_fcnPtr = &source_q_fcn_MDL1;
		source_u_fcnPtr = &source_u_fcn_MDL1;
		return;
	}
}

#endif /* I_CONDITIONS_H_ */
