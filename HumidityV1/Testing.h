/*
 * Testing.h
 *
 *  Created on: Oct 16, 2017
 *      Author: chuckjia
 */

#ifndef TESTING_H_
#define TESTING_H_
#include "TimeSteps.h"


/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Tests: Mesh Methods
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void printCellCenterCoord(int i, int j) {
	printf("Cell (%d, %d) center = (%f, %f)\n", i, j, getCellCenterX(i, j), getCellCenterP(i, j));
}

void printCellTopRightCoord(int i, int j) {
	printf("Cell (%d, %d) Top Right Vertex = (%f, %f)\n",
			i, j, getCellRightX(i, j), getCellTopRightP(i, j));
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Tests: Gaussian Elimination in the Projection Method
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Set test case 1: Nx = 6
void setTest1_gaussElimProj() {
	if (Nx != 6)
		throw "Error: Change the size of the mesh!";
	double test_a[5] = {2, 3, 4, 3, 3},
			test_b[5] = {3, 2, 3, 5, 4},
			test_c[5] = {2, 3, 1, 0, 3};
	for (int i = 1; i < Nx; ++i) {
		aInv_proj_cache[i] = 1 / test_a[i - 1];
		b_proj_cache[i] = test_b[i - 1];
		lambda_x_proj[i] = test_c[i - 1];
	}
}

// Set test case 2: Nx = 6
void setTest2_gaussElimProj() {
	if (Nx != 6)
		throw "Error: Change the size of the mesh!";
	double test_a[5] = {5, 6, 5, 6, 5},
			test_b[5] = {3, 3, 6, 9, 5},
			test_c[5] = {3, 2, 5, 1, 4};
	for (int i = 1; i < Nx; ++i) {
		aInv_proj_cache[i] = 1 / test_a[i - 1];
		b_proj_cache[i] = test_b[i - 1];
		lambda_x_proj[i] = test_c[i - 1];
	}
}

// Set test case 3: Nx = 10
void setTest3_gaussElimProj() {
	if (Nx != 10)
		throw "Error: Change the size of the mesh!";
	double test_a[9] = {5, 6, 5, 6, 5, 8, 9, 8, 6},
			test_b[9] = {3, 3, 6, 9, 5, 6, 7, 9, 6},
			test_c[9] = {3, 2, 5, 1, 3, 5, 6, 7, 8};
	for (int i = 1; i <= Nx; ++i) {
		aInv_proj_cache[i] = 1 / test_a[i - 1];
		b_proj_cache[i] = test_b[i - 1];
		lambda_x_proj[i] = test_c[i - 1];
	}
}

void selectTest_GaussElimProj() {
	setTest3_gaussElimProj();
}

void test_GaussElimProj() {
	try {
		selectTest_GaussElimProj();
	} catch (const char* msg) {
		cerr << msg << endl;
		return;
	}
	printf(">> Testing on the Gaussian Elimination in the Projection Method\n\n");

	// Print the original matrix
	printf("- The matrix to solve is:\n");
	for (int i = 1; i < Nx; ++i) {
		for (int j = 1; j <= Nx; ++j)
			if (i == j)
				printf("  %1.2f", 1 / aInv_proj_cache[i]);
			else if (i == j - 1)
				printf("  %1.2f", b_proj_cache[i]);
			else
				printf("  %1.2f", 0.);
		printf("  %1.2f\n", lambda_x_proj[i]);
	}
	for (int j = 1; j <= Nx; ++j)
		printf("  %1.2f", 1.);
	printf("  %1.2f\n", 0.);

	// Perform Gaussian elimination
	fillCache_d_proj();
	calcLambdax_gaussElim_proj();

	// Print result
	printf("\n- Result of the the Gaussian Elimination is\n[");
	for (int i = 1; i <= Nx; ++i)
		printf(" %1.15f;", lambda_x_proj[i]);
	printf(" ]\n");
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Tests: the Source Functions
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void testSourceFcns_helper(int numTestTimeSteps,
		double (*source_fcnPtr)(double, double, double, double, double, double, double)) {
	for (int tt = 0; tt < numTestTimeSteps; ++tt) {
		for (int i = 0; i <= Nx; ++i) {
			for (int j = 0; j < Np; ++j) {
				double T = T_sl[i][j], q = q_sl[i][j], u = u_sl[i][j], w = w_sl[i][j];
				double x = getCellCenterX(i, j), p = getCellCenterP(i, j), t = tt * Dt;
				printf("%1.2f ", (*source_T_fcnPtr)(T, q, u, w, x, p, t));
			}
			printf("\n");
		}
		printf("\n\n");
	}
}

void testSourceFcns() {
	int numTestTimeSteps = 5;
	printf("Testing on source function T:\n\n");
	testSourceFcns_helper(numTestTimeSteps, source_T_fcnPtr);
	printf("Testing on source function q:\n\n");
	testSourceFcns_helper(numTestTimeSteps, source_q_fcnPtr);
	printf("Testing on source function u:\n\n");
	testSourceFcns_helper(numTestTimeSteps, source_u_fcnPtr);
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Tests: the Initial Conditions
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void testIC() {
	enforceIC();
	peformAnalysis();
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Tests: the Quadrilateral Cell Interpolation Methods
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void printQuadCellCoefToFile() {
	FILE *f_a2 = fopen("Results/a2_quadCell.txt", "wb");
	FILE *f_a3 = fopen("Results/a3_quadCell.txt", "wb");
	FILE *f_a4 = fopen("Results/a4_quadCell.txt", "wb");
	for (int i = 0; i <= Nx; ++i)
		for (int j = 0; j <= Np; ++j) {
			fprintf(f_a2, "%1.20e ", a2_quadCell_cache[i][j]);
			fprintf(f_a3, "%1.20e ", a3_quadCell_cache[i][j]);
			fprintf(f_a4, "%1.20e ", a4_quadCell_cache[i][j]);
		}
	fclose(f_a2); fclose(f_a3); fclose(f_a4);
}

void printQuadCellDiagMatToFile() {
	FILE *f_e12 = fopen("Results/e12_diagMat_quadCell.txt", "wb");
	FILE *f_e22 = fopen("Results/e22_diagMat_quadCell.txt", "wb");
	for (int i = 0; i <= Nx; ++i)
		for (int j = 0; j <= Np; ++j) {
			fprintf(f_e12, "%1.20e ", e12_diagMatInv_quadCell_cache[i][j]);
			fprintf(f_e22, "%1.20e ", e22_diagMatInv_quadCell_cache[i][j]);
		}
	fclose(f_e12); fclose(f_e22);
	FILE *f = fopen("Results/e11e21_diagMat_quadCell.txt", "wb");
	fprintf(f, "%1.20e %1.20e ", e11_diagMatInv_quadCell_cache, e21_diagMatInv_quadCell_cache);
	fclose(f);
}

void testQuadCell() {
	writeParToFile();
	writeResToFile();
	printMeshToFile();
	printQuadCellCoefToFile();
	printQuadCellDiagMatToFile();
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Tests: Manufactured Solutions in Test Case 1
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void testExactSolnInMDL1_helper(double x, double p, double t) {
	double T_exactVal = exact_T_fcn_MDL1(x, p, t),
			T_xDer_exactVal = exact_T_xDer_fcn_MDL1(x, p, t),
			T_pDer_exactVal = exact_T_pDer_fcn_MDL1(x, p, t),
			T_tDer_exactVal = exact_T_tDer_fcn_MDL1(x, p, t),
			u_exactVal = exact_u_fcn_MDL1(x, p, t),
			u_xDer_exactVal = exact_u_xDer_fcn_MDL1(x, p, t),
			u_pDer_exactVal = exact_u_pDer_fcn_MDL1(x, p, t),
			u_tDer_exactVal = exact_u_tDer_fcn_MDL1(x, p, t),
			w_exactVal = exact_w_fcn_MDL1(x, p, t),
			w_pDer_exactVal = exact_w_pDer_fcn_MDL1(x, p, t);

	printf("\n>> Testing On the Manufactured Solutions in Test Case 1\n");
	printf("\nAt (x, p, t) = (%1.2f, %1.2f, %1.2f), the exact solutions are evaluated as\n",
			x, p, t);
	printf("\n- T = %1.5e, T_x = %1.5e, T_p = %1.5e, T_t = %1.5e\n",
			T_exactVal, T_xDer_exactVal, T_pDer_exactVal, T_tDer_exactVal);
	printf("- u = %1.5e, u_x = %1.5e, u_p = %1.5e, u_t = %1.5e\n",
			u_exactVal, u_xDer_exactVal, u_pDer_exactVal, u_tDer_exactVal);
	printf("- w = %1.5e, w_p = %1.5e\n", w_exactVal, w_pDer_exactVal);

	FILE *f = fopen("Results/test_exactSoln_MDL1.txt", "wb");
	fprintf(f,
			"%1.20e %1.20e %1.20e %1.20e %1.20e %1.20e %1.20e %1.20e %1.20e %1.20e %1.20e %1.20e %1.20e ",
			x, p, t,
			T_exactVal, T_xDer_exactVal, T_pDer_exactVal, T_tDer_exactVal,
			u_exactVal, u_xDer_exactVal, u_pDer_exactVal, u_tDer_exactVal,
			w_exactVal, w_pDer_exactVal);
	fclose(f);
}

void testExactSolnInMDL1() {
	double x = 34560.7273, p = 521.35, t = 132.35;
	testExactSolnInMDL1_helper(x, p, t);
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Tests: Upwind Method
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void print_r_uInterp_upwind_toFile() {
	FILE *f = fopen("Results/r_uInterp_upwind.txt", "wb");
	for (int i = 0; i <= Nx; ++i)
		fprintf(f, "%1.20e ", r_uInterp_cache[i]);
	fclose(f);
}

void emptyFcn() {
	// Empty
}

/*
 * Test A on upwind: fixed u and w
 */

void elimFluxes_u() {
	for (int i = 0; i <= Nx; ++i)
		for (int j = 0; j <= Np; ++j) {
			GG_u[i][j] = 0;
			FF_u[i][j] = 0;
		}
}

void calcFluxes_upwind_MDL3_test() {
	calcFluxes_upwind();
	elimFluxes_u();
}

void timeSteps_upwind_MDL3_test() {
	if (modelNo != 3)
		throw "Error: This test only works with Model 3!";
	calc_phix_fcnPtr = &emptyFcn;
	projU_fcnPtr = &emptyFcn;
	calc_w_fcnPtr = &emptyFcn;
	calcFluxes = &calcFluxes_upwind_MDL3_test;
	forwardEuler();
	//rk4();
}

/**
 * Test B on upwind: fixed u and w and reduced problem with no interpolation of u and w
 */

// Calculate all fluxes by upwind type Godunov scheme
void calcFluxes_upwind_reduced_MDL3_test() {
	// GG: fluxes on the TOP sides of cells
	for (int i = 1; i <= Nx; ++i)
		for (int j = 0; j <= Np; ++j) {
			double x = getCellCenterX(i, j), p = getCellTopLeftP(i, j);
			double wTopSideVal = exact_w_fcn_MDL3(x, p, 0);
			double velocity = Dx * wTopSideVal;
			if (velocity >= 0) {
				GG_T[i][j] =	 velocity * T_sl[i][j];
				GG_q[i][j] =	 velocity * q_sl[i][j];
			} else {
				GG_T[i][j] =	 velocity * T_sl[i][j + 1];
				GG_q[i][j] =	 velocity * q_sl[i][j + 1];
			}
		}

	// FF: fluxes on the RIGHT sides of cells
	for (int i = 0; i <= Nx; ++i)
		for (int j = 1; j <= Np; ++j) {
			double x = getCellRightX(i, j), p = getCellCenterP(i, j);
			double uRightSideVal = exact_u_fcn_MDL3(x, p, 0);
			double velocity = getCellRightDp(i, j) * uRightSideVal;
			if (velocity >= 0) {
				FF_T[i][j] =	 velocity * T_sl[i][j];
				FF_q[i][j] =	 velocity * q_sl[i][j];
			} else {
				FF_T[i][j] =	 velocity * T_sl[i + 1][j];
				FF_q[i][j] =	 velocity * q_sl[i + 1][j];
			}
		}
}

void timeSteps_upwind_reduced_MDL3_test() {
	if (modelNo != 3)
		throw "Error: This test only works with Model 3!";
	calc_phix_fcnPtr = &emptyFcn;
	projU_fcnPtr = &emptyFcn;
	calc_w_fcnPtr = &emptyFcn;
	calcFluxes = &calcFluxes_upwind_reduced_MDL3_test;
	//forwardEuler();
	rk4();
}

void testRun_MDL3() {
	try {
		//timeSteps_upwind_reduced_MDL3_test();
		timeSteps_upwind_MDL3_test();
	} catch (const char* msg) {
		cerr << msg << endl;
		return;
	}
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Tests on pB_xDer in the Physical Model
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// This test works with corresponding MATLAB script. To test, evaluate the function at the
// sample points in both c++ and MATLAB, and compare manually.
void test_pB_xDer_MDL0() {
	try {
		if (modelNo)
			throw "Error: This test only works for Model 0!";
	} catch (const char* msg) {
		cerr << msg << endl;
		return;
	}

	double arr[] = {1, 3e4, 3.5e4, 4e4, 4.5e4, 5e4, 6e4, 7e4};
	int n = sizeof(arr) / sizeof(arr[0]);
	for (int i = 0; i < n; ++i)
		printf("arr[%d] = %1.2e, pB_x = %1.7e \n", i, arr[i], pB_xDer_fcn_MDL0(arr[i]));
	printf("%f", c1_pBxDer_coef_MDL0);
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Tests: Decoupled System
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Enforce initial conditions on u for testing purposes
void enforceIC_exactVelocity(double t) {
	for (int i = 0; i < numCellsX; ++i) {
		double x = getCellCenterX(i);
		for (int j = 0; j < numCellsP; ++j) {
			double p = getCellCenterP(i, j);
			u_sl[i][j] = (*IC_u_fcnPtr)(x, p, t);
			w_sl[i][j] = (*IC_w_fcnPtr)(x, p, t);
		}
	}
	//enforceBC_topBD_numer_MDL1();
	//(*calc_w_fcnPtr)();
}

void rk4_decoupledVelocity_MDL102() {
	if (modelNo != 102)
		throw "This method only works with model 102";
	printf("\n- Running Runge-Kutta 4 method on time for decoupled velocity\n");
	int prog = -1;

	// The initial condition
	enforceIC();
	(*projU_fcnPtr)();
	(*calc_w_fcnPtr)();

	for (int tt = 0; tt < numTimeSteps; tt++) {
		// Print messages on calculation progress
		int progNew = tt * 100 / numTimeSteps;
		if (progNew > prog) {
			prog = progNew; printf("\r  - Current progress: %d%%", prog); fflush(stdout);
		}

		// Numerical calculation
		double t = Dt * tt;

		// RK4 Step 0
		copySoln(T_sl_copy, q_sl_copy, u_sl_copy);

		// RK4 Step 1
		enforceIC_exactVelocity(t);

		update_k_rk_fcnPtr = &update_k_rk_directAssign;
		preForwardEuler();
		forwardEuler_singleStep(t, halfDt, T_sl, q_sl, u_sl, ONE_SIXTH_CONST);
		postForwardEuler();

		// RK4 Step 2
		enforceIC_exactVelocity(t + halfDt);

		update_k_rk_fcnPtr = &update_k_rk_accum;
		preForwardEuler();
		forwardEuler_singleStep(t + halfDt, halfDt, T_sl_copy, q_sl_copy, u_sl_copy, ONE_THIRD_CONST);
		postForwardEuler();

		// RK4 Step 3
		enforceIC_exactVelocity(t + halfDt);

		update_k_rk_fcnPtr = &update_k_rk_accum;
		preForwardEuler();
		forwardEuler_singleStep(t + halfDt, Dt, T_sl_copy, q_sl_copy, u_sl_copy, ONE_THIRD_CONST);
		postForwardEuler();

		// RK4 Step 4
		enforceIC_exactVelocity(t + Dt);

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
	}
	printf("\r  - Runge-Kutta 4 method complete\n");
}

void testRun_MDL102() {
	try {
		rk4_decoupledVelocity_MDL102();
	} catch (const char* msg) {
		cerr << msg << endl;
		exit(EXIT_FAILURE);
	}
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Use Exact w but Numerically compute u
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void calc_w_exact(double t) {
	for (int i = 0; i < numCellsX; ++i) {
		double x = getCellCenterX(i);
		for (int j = 0; j < Np; ++j) {
			double p = getCellCenterP(i, j);
			w_sl[i][j] = (*IC_w_fcnPtr)(x, p, t);
		}
	}
	enforceNonPenetrationBC_topBD_numer();
}

// Use exact w but numerically enforce non-penetration BC on the top for w
void rk4_exactW_withTopBC_MDL1_test() {
	if (modelNo != 1)
		throw "Error: Only works for Test Case 1!";
	printf("\n- Running Runge-Kutta 4 method on time with exact w\n");
	int prog = -1;

	// The initial condition
	enforceIC();
	(*projU_fcnPtr)();

	for (int tt = 0; tt < numTimeSteps; tt++) {
		// Print messages on calculation progress
		int progNew = tt * 100 / numTimeSteps;
		if (progNew > prog) {
			prog = progNew; printf("\r  - Current progress: %d%%", prog); fflush(stdout);
		}

		// Numerical calculation
		double t = Dt * tt;

		// RK4 Step 0
		copySoln(T_sl_copy, q_sl_copy, u_sl_copy);

		// RK4 Step 1
		calc_w_exact(t);

		update_k_rk_fcnPtr = &update_k_rk_directAssign;
		preForwardEuler();
		forwardEuler_singleStep(t, halfDt, T_sl, q_sl, u_sl, ONE_SIXTH_CONST);
		postForwardEuler();


		// RK4 Step 2
		calc_w_exact(t + halfDt);

		update_k_rk_fcnPtr = &update_k_rk_accum;
		preForwardEuler();
		forwardEuler_singleStep(t + halfDt, halfDt, T_sl_copy, q_sl_copy, u_sl_copy, ONE_THIRD_CONST);
		postForwardEuler();

		// RK4 Step 3
		calc_w_exact(t + halfDt);

		update_k_rk_fcnPtr = &update_k_rk_accum;
		preForwardEuler();
		forwardEuler_singleStep(t + halfDt, Dt, T_sl_copy, q_sl_copy, u_sl_copy, ONE_THIRD_CONST);
		postForwardEuler();

		// RK4 Step 4
		calc_w_exact(t + Dt);

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
	}
	printf("\r  - Runge-Kutta 4 method complete\n");
}

void testRun_exactW_MDL1() {
	try {
		rk4_exactW_withTopBC_MDL1_test();
	} catch (const char* msg) {
		cerr << msg << endl;
		exit(EXIT_FAILURE);
	}
}

void printPhixErrToFile(double t) {
	FILE *f = fopen("Results/phix_err.txt", "wb");
	for (int i = 1; i <= Nx; ++i) {
		double x = getCellCenterX(i), DpVal = getCellCenterDp(i), sum = 0;
		for (int j = 1; j <= Np; ++j) {
			double p = getCellCenterP(i, j), pInCurrCell = p - (j - 1) * DpVal;
			double T_xDer_val = exact_T_xDer_fcn_MDL1(x, p, t);
			double phixVal = sum - R_CONST * T_xDer_val / p * pInCurrCell;
			sum -= R_CONST * T_xDer_val / p * DpVal;
			fprintf(f, "%1.20e ", phixVal - phix_sl[i][j]);
		}
	}
	fclose(f);
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Testing
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void testing() {
	runTimeSteps();
	//testRun_exactW_MDL1();
	// testRun_MDL102();
	peformAnalysis();
	printPhixErrToFile(finalTime);
}
#endif /* TESTING_H_ */
