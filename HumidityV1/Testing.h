/*
 * Testing.h
 *
 *  Created on: Oct 16, 2017
 *      Author: chuckjia
 */

#ifndef TESTING_H_
#define TESTING_H_
#include <iostream>
#include "Analysis.h"
using namespace std;

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Testing Mesh Methods
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void printMeshToFile() {
	FILE *fGridX = fopen("gridX.txt", "wb");
	FILE *fGridP = fopen("gridP.txt", "wb");
	for (int i = 0; i < numGridPtsX; ++i)
		for (int j = 0; j < numGridPtsP; ++j) {
			fprintf(fGridX, "%f ", getCellLeftX(i, j));
			fprintf(fGridP, "%f ", getCellBottLeftP(i, j));
		}
	fclose(fGridX);
	fclose(fGridP);

	FILE *fCellCenterX = fopen("cellCenterX.txt", "wb");
	FILE *fCellCenterP = fopen("cellCenterP.txt", "wb");
	FILE *fCellVol = fopen("cellVol.txt", "wb");
	for (int i = 0; i < numCellsX; ++i)
		for (int j = 0; j < numCellsP; ++j) {
			fprintf(fCellCenterX, "%f ", getCellCenterX(i, j));
			fprintf(fCellCenterP, "%f ", getCellCenterP(i, j));
			fprintf(fCellVol, "%f ", getCellVol(i, j));
		}
	fclose(fCellCenterX);
	fclose(fCellCenterP);
	fclose(fCellVol);
}

void printMeshPar() {
	FILE *f = fopen("par.txt", "wb");
	fprintf(f, "%f %f %f %f %f %d %d",
			x0, xf, pA, (*pB_fcnPtr)(x0), (*pB_fcnPtr)(xf), Nx, Np);
	fclose(f);
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
	gaussElim_proj();

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
 * Testing
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void testing() {
	timeSteps();
	peformAnalysis();
}
#endif /* TESTING_H_ */
