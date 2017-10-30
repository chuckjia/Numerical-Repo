/*
 * Testing.h
 *
 *  Created on: Oct 16, 2017
 *      Author: chuckjia
 */

#ifndef TESTING_H_
#define TESTING_H_
#include <iostream>
#include "TimeSteps.h"
using namespace std;

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Testing Mesh Methods
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void printCellCenterCoord(int i, int j) {
	printf("Cell (%d, %d) center = (%f, %f)\n", i, j, getCellCenterX(i, j), getCellCenterP(i, j));
}

void printCellTopRightCoord(int i, int j) {
	printf("Cell (%d, %d) Top Right Vertex = (%f, %f)\n",
			i, j, getCellRightX(i, j), getCellTopRightP(i, j));
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Testing the Gaussian Elimination in the Projection Method
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
 * Testing the Source Functions
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void testSourceFcns_helper(int numTestTimeSteps,
		double (*source_fcnPtr)(double, double, double, double, double, double)) {
	for (int tt = 0; tt < numTestTimeSteps; ++tt) {
		for (int i = 0; i <= Nx; ++i) {
			for (int j = 0; j < Np; ++j) {
				double T = T_sl[i][j], q = q_sl[i][j], u = u_sl[i][j];
				double x = getCellCenterX(i, j), p = getCellCenterP(i, j), t = tt * Dt;
				printf("%1.2f ", (*source_T_fcnPtr)(T, q, u, x, p, t));
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
 * Testing the Initial Conditions
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void testIC() {
	enforceIC();
	peformAnalysis();
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Testing the Quadrilateral Cell Interpolation Methods
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
	printParToFile();
	writeResToFile();
	printMeshToFile();
	printQuadCellCoefToFile();
	printQuadCellDiagMatToFile();
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Testing On the Manufactured Solutions in Test Case 1
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

void writeExactSolnToFile() {
	// Write final numerical solution to files
	FILE *exact_T = fopen("Results/T_exact.txt", "wb");
	FILE *exact_q = fopen("Results/q_exact.txt", "wb");
	FILE *exact_u = fopen("Results/u_exact.txt", "wb");
	FILE *exact_w = fopen("Results/w_exact.txt", "wb");
	double t = finalTime;
	for (int i = 0; i < numCellsX; ++i) {
		double x = getCellCenterX(i);
		for (int j = 0; j < numCellsP; ++j) {
			double p = getCellCenterP(i, j);
			fprintf(exact_T, "%f ", exact_T_fcn_MDL1(x, p, t));
			fprintf(exact_q, "%f ", exact_q_fcn_MDL1(x, p, t));
			fprintf(exact_u, "%f ", exact_u_fcn_MDL1(x, p, t));
			fprintf(exact_w, "%f ", exact_w_fcn_MDL1(x, p, t));
		}
	}
	fclose(exact_T); fclose(exact_q); fclose(exact_u); fclose(exact_w);
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Testing
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void testing() {
	timeSteps();
	peformAnalysis();
	writeExactSolnToFile();
	// testExactSolnInMDL1();
}
#endif /* TESTING_H_ */
