/*
 * Analysis.h
 *
 *  Created on: Oct 26, 2017
 *      Author: chuckjia
 */

#ifndef ANALYSIS_H_
#define ANALYSIS_H_
#include "Conditions.h"

const double GRTD_PREC_CONST = 1e-15;  // Guaranteed precision

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Print Messages
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void printTitle() {
	printf("===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== \n");
	printf(">> Solving Atmospherical Model Using Upwind-type Godunov Scheme\n");
	printf("===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== \n\n");
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Print All Parameters For Diagnostics
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void printDiagnostics() {
	printf("\n>> DIAGNOSTIC INFORMATION\n");
	printf("\nParameters: \n");

	printf("\n - Time\n");
	printf("    [1] Dt = %1.2e, finalTime = %1.2es, numTimeSteps = %d\n", Dt, finalTime, numTimeSteps);

	printf("\n - Domain geometry\n");
	printf("    [1] x0 = %1.2f, xf = %1.2f \n", x0, xf);
	printf("    [2] pA = %1.2f, pB(x0) = %1.2f, pB[(x0+xf)/2] = %1.2f, pB(xf) = %1.2f\n",
			pA, (*pB_fcnPtr)(x0), (*pB_fcnPtr)(0.5 * (x0 + xf)), (*pB_fcnPtr)(xf));

	printf("\n - Mesh specifications\n");
	printf("    [1] Nx = %d, Np = %d\n", Nx, Np);
	printf("    [2] Dx = %1.2f\n", Dx);
	printf("    [3] Dp(x0) = %1.2f, Dp[(x0+xf)/2] = %1.2f, Dp(xf) = %1.2f\n",
			getCellLeftDp(1), getCellLeftDp(lastRealIndexX / 2), getCellLeftDp(lastRealIndexX));
}

void printSchemeDescription() {
	printf("\n>> SCHEME DESCRIPTION\n");
	printf("\n  MODEL %d:  \n", modelNo);

	printf("  [1] Time:  ");
	printf("Dt = %1.2e, finalTime = %1.2es, numTimeSteps = %d\n", Dt, finalTime, numTimeSteps);

	printf("  [2] Mesh Specs:  ");
	printf("Nx = %d, Np = %d\n", Nx, Np);
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Print The Mesh To File
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void printMeshToFile() {
	FILE *fGridX = fopen("Results/meshGridX.txt", "wb");
	FILE *fGridP = fopen("Results/meshGridP.txt", "wb");
	for (int i = 0; i < numGridPtsX; ++i)
		for (int j = 0; j < numGridPtsP; ++j) {
			fprintf(fGridX, "%1.20e ", getCellLeftX(i, j));
			fprintf(fGridP, "%1.20e ", getCellBottLeftP(i, j));
		}
	fclose(fGridX);
	fclose(fGridP);

	FILE *fCellCenterX = fopen("Results/cellCentersX.txt", "wb");
	FILE *fCellCenterP = fopen("Results/cellCentersP.txt", "wb");
	FILE *fCellVol = fopen("Results/cellVol.txt", "wb");
	FILE *fCellSideLen = fopen("Results/cellSideLen.txt", "wb");
	for (int i = 0; i < numCellsX; ++i)
		for (int j = 0; j < numCellsP; ++j) {
			fprintf(fCellCenterX, "%1.20e ", getCellCenterX(i, j));
			fprintf(fCellCenterP, "%1.20e ", getCellCenterP(i, j));
			fprintf(fCellVol, "%1.20e ", getCellVol(i, j));
			fprintf(fCellSideLen, "%1.20e ", getCellBottSideLen(i, j));
		}
	fclose(fCellCenterX);
	fclose(fCellCenterP);
	fclose(fCellVol);
	fclose(fCellSideLen);

	FILE *fCellTopNormX = fopen("Results/cellTopNormX.txt", "wb");
	FILE *fCellTopNormP = fopen("Results/cellTopNormP.txt", "wb");
	for (int i = 1; i <= Nx; ++i)
		for (int j = 0; j <= Np; ++j) {
			fprintf(fCellTopNormX, "%1.20e ", getCellTopSideNormX(i, j));
			fprintf(fCellTopNormP, "%1.20e ", getCellTopSideNormP(i, j));
		}
	fclose(fCellTopNormX);
	fclose(fCellTopNormP);
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Calculate and Print Errors
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

double cellErr_helper(double numericalVal, double exactVal, double vol) {
	return vol * pow(exactVal - numericalVal, 2);
}

double relativeL2Err_helper(double num, double denom) {
	if (denom > GRTD_PREC_CONST)
		return sqrt(num / denom);
	else if (num < GRTD_PREC_CONST)
		return 0;
	return -99;
}

// Print the L2 relative and absolute errors
void showL2Errors(double t) {
	double num_T = 0, num_q = 0, num_u = 0, num_w = 0,
			denom_T = 0, denom_q = 0, denom_u = 0, denom_w = 0;
	for (int i = 1; i < lastGhostIndexX; ++i) {
		double x = getCellCenterX(i);
		for (int j = 1; j < lastGhostIndexP; ++j) {
			double p = getCellCenterP(i, j), vol = getCellVol(i, j);
			double exactVal_T = (*IC_T_fcnPtr)(x, p, t),
					exactVal_q = (*IC_q_fcnPtr)(x, p, t),
					exactVal_u = (*IC_u_fcnPtr)(x, p, t),
					exactVal_w = (*IC_w_fcnPtr)(x, p, t);
			// Calculate errors
			num_T += cellErr_helper(exactVal_T, T_sl[i][j], vol);
			num_q += cellErr_helper(exactVal_q, q_sl[i][j], vol);
			num_u += cellErr_helper(exactVal_u, u_sl[i][j], vol);
			num_w += cellErr_helper(exactVal_w, w_sl[i][j], vol);
			denom_T += vol * pow(exactVal_T, 2);
			denom_q += vol * pow(exactVal_q, 2);
			denom_u += vol * pow(exactVal_u, 2);
			denom_w += vol * pow(exactVal_w, 2);
		}
	}
	// double absL2Err_T = sqrt(num_T), absL2Err_q = sqrt(num_q),
	// absL2Err_u = sqrt(num_u), absL2Err_w(num_w);
	double relativeL2Err_T = relativeL2Err_helper(num_T, denom_T),
			relativeL2Err_q = relativeL2Err_helper(num_q, denom_q),
			relativeL2Err_u = relativeL2Err_helper(num_u, denom_u),
			relativeL2Err_w = relativeL2Err_helper(num_w, denom_w);

	printf("\n  - L2 relative error for T, q, u, w = [%1.4e, %1.4e, %1.4e, %1.4e]\n",
			relativeL2Err_T, relativeL2Err_q, relativeL2Err_u, relativeL2Err_w);
	printf("  - L2 discrete norm for T, q, u, w = [%1.4e, %1.4e, %1.4e, %1.4e]\n",
			sqrt(denom_T), sqrt(denom_q), sqrt(denom_u), sqrt(denom_w));
	// printf("\n  - L2 absolute error for T, q, u, w = [%1.4e, %1.4e, %1.4e, %1.4e]\n",
	// absL2Err_T, absL2Err_q, absL2Err_u, absL2Err_w);
}

void showL2Errors() {
	printf("\n>> NUMERICAL ERRORS\n");
	showL2Errors(finalTime);
}

void writeParToFile() {
	FILE *f = fopen("Results/par.txt", "wb");
	fprintf(f, "%1.15f %1.15f %1.15f %1.15f %1.15f %d %d %1.15f %d",
			x0, xf, pA, (*pB_fcnPtr)(x0), (*pB_fcnPtr)(xf), Nx, Np, Dt, numTimeSteps);
	fclose(f);
}

// Write the result to file: res.txt for the solution and err.txt for the error
// This function is used in plotting the solution and the error
void writeResToFile() {
	// Write final numerical solution to files
	FILE *res_T = fopen("Results/T_soln.txt", "wb");
	FILE *res_q = fopen("Results/q_soln.txt", "wb");
	FILE *res_u = fopen("Results/u_soln.txt", "wb");
	FILE *res_w = fopen("Results/w_soln.txt", "wb");
	for (int i = 0; i < numCellsX; ++i)
		for (int j = 0; j < numCellsP; ++j) {
			fprintf(res_T, "%f ", T_sl[i][j]);
			fprintf(res_q, "%f ", q_sl[i][j]);
			fprintf(res_u, "%f ", u_sl[i][j]);
			fprintf(res_w, "%f ", w_sl[i][j]);
		}
	fclose(res_T); fclose(res_q); fclose(res_u); fclose(res_w);

	// Write final numerical errors to files
	FILE *err_T = fopen("Results/T_err.txt", "wb");
	FILE *err_q = fopen("Results/q_err.txt", "wb");
	FILE *err_u = fopen("Results/u_err.txt", "wb");
	FILE *err_w = fopen("Results/w_err.txt", "wb");
	for (int i = 0; i < numCellsX; ++i) {
		double x = getCellCenterX(i);
		for (int j = 0; j < numCellsP; ++j) {
			double p = getCellCenterP(i, j);
			double exactSoln_T = (*IC_T_fcnPtr)(x, p, finalTime),
					exactSoln_q = (*IC_q_fcnPtr)(x, p, finalTime),
					exactSoln_u = (*IC_u_fcnPtr)(x, p, finalTime),
					exactSoln_w = (*IC_w_fcnPtr)(x, p, finalTime);
			fprintf(err_T, "%f ", T_sl[i][j] - exactSoln_T);
			fprintf(err_q, "%f ", q_sl[i][j] - exactSoln_q);
			fprintf(err_u, "%f ", u_sl[i][j] - exactSoln_u);
			fprintf(err_w, "%f ", w_sl[i][j] - exactSoln_w);
		}
	}
	fclose(err_T); fclose(err_q); fclose(err_u); fclose(err_w);
}

void peformAnalysis() {
	clock_t start = clock();

	printSchemeDescription();
	showL2Errors();
	writeParToFile();
	printMeshToFile();
	writeResToFile();

	printf("\n- Analysis complete. Time used = %1.2fs.\n",
			((double) (clock() - start)) / CLOCKS_PER_SEC);
}

#endif /* ANALYSIS_H_ */
