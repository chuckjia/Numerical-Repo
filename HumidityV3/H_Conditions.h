/*
 * Conditions.h
 *
 *  Created on: Oct 14, 2017
 *      Author: chuckjia
 *  !!AlphaVersion!!
 */

#ifndef H_CONDITIONS_H_
#define H_CONDITIONS_H_
#include "G_WPhix.h"

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Boundary Conditions
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Function pointer: function to enforce all boundary conditions in one model
void (*enforceBC_fptr)();

/* ----- ----- ----- ----- ----- -----
 * Common boundary conditions
 * ----- ----- ----- ----- ----- ----- */

/*
 * In implementing all the BCs, we assume pB is sufficiently flat on the left and right domain
 * sides. This guarantees the left and right ghost cells are not distorted.
 */

/**
 * LEFT side of domain
 */

// Enforce Dirichlet BCs on the LEFT boundary: with boundary values specified by a function, defined by (*bdVal_fptr)(p)
void enforceDirichlet_leftBD(double sl[numCellX][numCellP], double (&bdVal_fptr)(double)) {
	for (int j = 1; j <= Np; ++j) {
		// This interpolation on p only works for pB() that are sufficiently flat on the left side of domain
		double p = 0.5 * (getCellBottLeftP(1, j) + getCellTopLeftP(1, j));
		sl[0][j] = 2 * bdVal_fptr(p) - sl[1][j];
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
		sl[0][j] = -sl[1][j];
}

// Cache to store the boundary values specified on the left side of boundary
double T_leftBdVal_[numCellP], q_leftBdVal_[numCellP], u_leftBdVal_[numCellP];

// Enforce Dirichlet BCs on the LEFT boundary: with boundary value from cache
void enforceDirichlet_leftBD(double sl[numCellX][numCellP], double bdVal_[numCellP]) {
	for (int j = 1; j <= Np; ++j)
		sl[0][j] = 2 * bdVal_[j] - sl[1][j];
	// sl[0][j] = bdVal_[j];  // CHANGED!
}

// Enforce Neumann BCs on the left boundary
void enforceNeumann_leftBD(double sl[numCellX][numCellP]) {
	for (int j = 1; j <= Np; ++j)
		sl[0][j] = sl[1][j];
}

/**
 * RIGHT side of domain
 */

// Enforce Dirichlet BCs on the RIGHT boundary: with boundary value 0
void enforceDirichlet_rightBD(double sl[numCellX][numCellP]) {
	for (int j = 1; j <= Np; ++j)
		sl[lastGhostIndexX][j] = -sl[Nx][j];
}

// Enforce Neumann BCs on the right boundary
void enforceNeumann_rightBD(double sl[numCellX][numCellP]) {
	for (int j = 1; j <= Np; ++j)
		sl[lastGhostIndexX][j] = sl[Nx][j];
}

/**
 * BOTTOM side of domain
 */

// Enforce Dirichlet BCs on the BOTTOM boundary: with boundary value 0
void enforceDirichlet_bottBD(double sl[numCellX][numCellP]) {
	for (int i = 1; i <= Nx; ++i)
		sl[i][0] = -sl[i][1];
}

/**
 * TOP side of domain
 */

// !!AlphaVersion!!
// Enforce Non-penetration boundary condition on the top side of the domain. Enforcement is done by extrapolation mathematically
void enforceNonPenetrationBC_topBD_math() {
	for (int i = 1; i <= Nx; ++i) {
		double x = getCellCenterX(i);
		w_[i][lastGhostIndexP] = (*pBxDer_fptr)(x) * (u_[i][Np] + u_[i][lastGhostIndexP]) - w_[i][Np];
	}
}

// !!AlphaVersion!!
// Enforce Non-penetration boundary condition on the top side of the domain
void enforceNonPenetrationBC_topBD() {
	//	// This method ensures the flux only comes from within the domain, not from the ghost cells, by mirroring within domain boundary cell
	//	// values to ghost cells
	//	for (int i = 1; i <= Nx; ++i) {
	//		T_[i][lastGhostIndexP] = T_[i][Np];
	//		q_[i][lastGhostIndexP] = q_[i][Np];
	//		u_[i][lastGhostIndexP] = u_[i][Np];
	//		// w_[i][lastGhostIndexP] = w_[i][Np];
	//	}

	// This part ensures the upwind scheme only draws flux values from within the domain, not from ghost cells
	for (int i = 1; i <= Nx; ++i) {
		double norm_x = getCellTopSideNormVecX(i, Np), norm_p = getCellTopSideNormVecP(i, Np);
		w_[i][lastGhostIndexP] = -norm_x * (u_[i][Np] + u_[i][lastGhostIndexP]) / norm_p - w_[i][Np];
	}
}

/* ----- ----- ----- ----- ----- -----
 * Physical Case 0: BCs
 * ----- ----- ----- ----- ----- ----- */

// Helpers for projecting u values on the left side of domain at initial time t = 0
double helper1_fillCache_leftBdVal_MDL0(double p) { return sin(M_PI * p / p0_CONST); }
double helper2_fillCache_leftBdVal_MDL0(double x) { return cos(TWO_PI * _n_initU_MDL0 / xf * x); }
void enforceIC();  // Initialize function as it will be used in fillCache_leftBdVal_MDL0()

// !!AlphaVersion!!
// Pre-calculate the left boundary values and store in cache
void fillCache_leftBdVal_MDL0() {
	/*// The following is to project u for the left boundary values
	enforceIC();
	// Calculate lambda_x(x0)
	double pB = pB_fcn_MDL0(x0), x1 = getCellCenterX(1),
			lambda_x_leftBdVal =  lambdax_proj_[1] -
			2 * p0_CONST / (M_PI * (pB - pA))
	 * (helper1_fillCache_leftBdVal_MDL0(pB) - helper1_fillCache_leftBdVal_MDL0(pA))
	 * (helper2_fillCache_leftBdVal_MDL0(x1) - helper2_fillCache_leftBdVal_MDL0(x0));*/

	for (int j = 0; j < numCellP; ++j) {
		double p = getCellCenterP(1, j), T = init_T_fcn_MDL0(x0, p, 0);
		T_leftBdVal_[j] = T;  // Put T values in cache
		q_leftBdVal_[j] = qs_fcn(T, p) * 0.5;  // Put q values in cache
		u_leftBdVal_[j] = init_u_fcn_MDL0(0, p, 0); // - lambda_x_leftBdVal;  // Put u values in cache
	}
}

// Enforce all boundary conditions
void enforceBC_MDL0() {
	// Left boundary: Dirichlet BC
	enforceDirichlet_leftBD(T_, T_leftBdVal_);
	enforceDirichlet_leftBD(q_, q_leftBdVal_);
	enforceDirichlet_leftBD(u_, u_leftBdVal_);
	enforceNeumann_leftBD(w_);

	// Right boundary: Neumann BC
	enforceNeumann_rightBD(T_);
	enforceNeumann_rightBD(q_);
	enforceNeumann_rightBD(u_);
	enforceNeumann_rightBD(w_);

	// Bottom boundary: for w and phi_x, Dirichlet BC with boudnary value 0
	// Conditions are enforced directly in the calculation of the two solutions

	// Top boundary: for u and w
	if (_enforceTopBC_)
		enforceNonPenetrationBC_topBD();  // CHANGED!
	// enforceNonPenetrationBC_topBD_math();
}

/* ----- ----- ----- ----- ----- -----
 * Test Case 1 BCs
 * ----- ----- ----- ----- ----- ----- */

double leftBdVal_T_fcn_MDL1(double p) { return exact_T_fcn_MDL1(0, p, 0); }

void enforceBC_MDL1() {
	// Left boundary
	enforceDirichlet_leftBD(T_, leftBdVal_T_fcn_MDL1);
	//	enforceDirichlet_leftBD(T_);
	enforceDirichlet_leftBD(q_);
	enforceDirichlet_leftBD(u_);

	// Right boundary: Neumann BC
	enforceNeumann_rightBD(T_);
	enforceNeumann_rightBD(q_);
	enforceNeumann_rightBD(u_);

	// Bottom boundary: for w, Dirichlet BC with boudnary value 0

	// Top boundary: for u and w
	//	enforceNonPenetrationBC_topBD();
}


/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Source Functions
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Function pointers: mathematical functions to calculate source term values
double (*source_T_fptr)(double T, double q, double u, double w, double x, double p, double t),
		(*source_q_fptr)(double T, double q, double u, double w, double x, double p, double t),
		(*source_u_fptr)(double T, double q, double u, double w, double x, double p, double t);


/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Initial Conditions
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Function pointers: mathematical functions to calculate initial conditions
// We implement input t for testing purposes. The initial condition for w is only for testing purposes
double (*initT_fptr)(double x, double p, double t), (*initQ_fptr)(double x, double p, double t),
		(*initU_fptr)(double x, double p, double t), (*initW_fptr)(double x, double p, double t);

// Enforce initial conditions at cell centers
void enforceIC() {
	// We do NOT set initial conditions for w, as in the model
	for (int i = 0; i < numCellX; ++i)
		for (int j = 0; j < numCellP; ++j) {
			double x = getCellCenterX(i, j), p = getCellCenterP(i, j);
			T_[i][j] = (*initT_fptr)(x, p, 0);
			q_[i][j] = (*initQ_fptr)(x, p, 0);
			u_[i][j] = (*initU_fptr)(x, p, 0);
			// w_[i][j] = (*initW_fptr)(x, p, 0);  // Changed
		}
	(*projU_fptr)();
	(*calcW_fptr)();
}


void enforceIC(double sl[numCellX][numCellP], double (*initSl_fptr)(double, double, double), double t) {
	for (int i = 0; i < numCellX; ++i)
		for (int j = 0; j < numCellP; ++j)
			sl[i][j] = (*initSl_fptr)(getCellCenterX(i, j), getCellCenterP(i, j), t);
}


/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Averaging Method
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void (*aveSoln_fptr)(int tt);  // Function pointer to averaging method
bool aveMethodApplied = true;

// Average only a single solution
// !!AlphaVersion!! This implementation uses information in the ghost cells
void aveSoln(double sl[numCellX][numCellP]) {
	for (int i = 1; i <= Nx; ++i)
		for (int j = 1; j <= Np; ++j)
			sl[i][j] = 0.5 * (sl[i][j] + sl[i - 1][j]);
}

void aveSoln_vertical(double sl[numCellX][numCellP]) {
	for (int i = 1; i <= Nx; ++i)
		for (int j = 2; j <= Np; ++j)
			sl[i][j] = 0.5 * (sl[i][j - 1] + sl[i][j]);
}

// Apply average method on all solutions
void aveSoln(int tt) {
	if (tt == 0) return;

	if (!(tt % aveSolnFreq_T))     aveSoln(T_);
	if (!(tt % aveSolnFreq_q))     aveSoln(q_);
	if (!(tt % aveSolnFreq_u))     aveSoln(u_);
	if (!(tt % aveSolnFreq_w))     aveSoln(w_);
	if (!(tt % aveSolnFreq_phix))  aveSoln(phix_);

	(*enforceBC_fptr)();
}


/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Set All Conditions: ICs, BCs and Source Functions
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void setConditions() {
	switch (modelNo) {
	case 0:
		// Initial conditions
		initT_fptr = &init_T_fcn_MDL0;
		initQ_fptr = &init_q_fcn_MDL0;
		initU_fptr = &init_u_fcn_MDL0;
		initW_fptr = &zero_fcn;
		// Boundary conditions
		fillCache_leftBdVal_MDL0();
		enforceBC_fptr = &enforceBC_MDL0;
		// Source functions
		source_T_fptr = &source_T_fcn_MDL0;
		source_q_fptr = &source_q_fcn_MDL0;
		source_u_fptr = &source_u_fcn_MDL0;
		// Averaging method
		aveSoln_fptr = &aveSoln;
		aveMethodApplied = true;
		break;

	case 1:
		// Initial conditions
		initT_fptr = &exact_T_fcn_MDL1;
		initQ_fptr = &exact_q_fcn_MDL1;
		initU_fptr = &exact_U_fcn_MDL1;
		initW_fptr = &exact_W_fcn_MDL1;
		// Boundary conditions
		enforceBC_fptr = &enforceBC_MDL1;
		// Source functions
		source_T_fptr = &source_T_fcn_MDL1;
		source_q_fptr = &source_q_fcn_MDL1;
		source_u_fptr = &source_u_fcn_MDL1;
		// Averaging method
		aveSoln_fptr = &empty_fcn;
		aveMethodApplied = false;
		break;
	}
}

#endif /* H_CONDITIONS_H_ */
