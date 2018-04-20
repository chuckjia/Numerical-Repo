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
double computationTime = 0;
double relative_L2Err_T_global, relative_L2Err_q_global,
relative_L2Err_u_global, relative_L2Err_w_global;
FILE *T_norm_filePtr, *q_norm_filePtr, *u_norm_filePtr, *w_norm_filePtr, *time_filePtr;

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
	printf("    [1] Dt = %1.2e, finalTime = %1.2es, numTimeSteps = %d\n", Dt, finalTime, numTimeStep);

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
	printf("Dt = %1.2e, finalTime = %1.2es, numTimeSteps = %d\n", Dt, finalTime, numTimeStep);

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
	return -99999e10;
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
	double absL2Err_T = sqrt(num_T), absL2Err_q = sqrt(num_q),
			absL2Err_u = sqrt(num_u), absL2Err_w(num_w);
	double relativeL2Err_T = relativeL2Err_helper(num_T, denom_T),
			relativeL2Err_q = relativeL2Err_helper(num_q, denom_q),
			relativeL2Err_u = relativeL2Err_helper(num_u, denom_u),
			relativeL2Err_w = relativeL2Err_helper(num_w, denom_w);

	printf("\n  - L2 relative error for T, q, u, w = [%1.7e, %1.7e, %1.7e, %1.7e]\n",
			relativeL2Err_T, relativeL2Err_q, relativeL2Err_u, relativeL2Err_w);
	printf("  - L2 discrete norm for T, q, u, w = [%1.7e, %1.7e, %1.7e, %1.7e]\n",
			sqrt(denom_T), sqrt(denom_q), sqrt(denom_u), sqrt(denom_w));
	printf("  - L2 absolute error for T, q, u, w = [%1.7e, %1.7e, %1.7e, %1.7e]\n",
			absL2Err_T, absL2Err_q, absL2Err_u, absL2Err_w);
	relative_L2Err_T_global = relativeL2Err_T;
	relative_L2Err_q_global = relativeL2Err_q;
	relative_L2Err_u_global = relativeL2Err_u;
	relative_L2Err_w_global = relativeL2Err_w;
}

void showL2Errors() {
	printf("\n>> NUMERICAL ERRORS\n");
	showL2Errors(finalTime);
}

// Print the L2 relative and absolute errors
double calcL2Norm(double sl[numCellsX][numCellsP], double t) {
	double norm = 0;
	for (int i = 1; i <= Nx; ++i) {
		double vol = getCellVol(i);
		for (int j = 1; j <= Np; ++j)
			norm += vol * pow(sl[i][j], 2);
	}
	return sqrt(norm);
}

void calcL2Norm(double t) {
	double norm_T = 0, norm_q = 0, norm_u = 0, norm_w = 0;
	for (int i = 1; i <= Nx; ++i) {
		double vol = getCellVol(i);
		for (int j = 1; j <= Np; ++j) {
			norm_T += vol * pow(T_sl[i][j], 2);
			norm_q += vol * pow(q_sl[i][j], 2);
			norm_u += vol * pow(u_sl[i][j], 2);
			norm_w += vol * pow(w_sl[i][j], 2);
		}
	}
	norm_T = sqrt(norm_T); norm_q = sqrt(norm_q); norm_u = sqrt(norm_u); norm_w = sqrt(norm_w);
	printf("\n    L2 norm (T, q, u, w) = (%1.2f, %1.2f, %1.2f, %1.2f)",
			norm_T - 2.055e6, norm_q, norm_u, norm_w);
	fprintf(T_norm_filePtr, "%1.20e ", norm_T);
	fprintf(q_norm_filePtr, "%1.20e ", norm_q);
	fprintf(u_norm_filePtr, "%1.20e ", norm_u);
	fprintf(w_norm_filePtr, "%1.20e ", norm_w);
	fprintf(time_filePtr, "%1.20e ", t);
}

void printParamToFile() {
	FILE *f = fopen("Results/par.txt", "wb");
	fprintf(f, "%1.20e %1.20e %1.20e %1.20e %1.20e %d %d %1.20e %d %d %d",
			x0, xf, pA, (*pB_fcnPtr)(x0), (*pB_fcnPtr)(xf), Nx, Np, Dt, numTimeStep,
			movieFrameFreq, aveFreq);
	fclose(f);
}

// Write the result to file: res.txt for the solution and err.txt for the error
// This function is used in plotting the solution and the error
void printResToFile() {
	// Numerical solution files
	FILE *res_T = fopen("Results/T_soln.txt", "wb");
	FILE *res_q = fopen("Results/q_soln.txt", "wb");
	FILE *res_u = fopen("Results/u_soln.txt", "wb");
	FILE *res_w = fopen("Results/w_soln.txt", "wb");
	// Numerical errors files
	FILE *err_T = fopen("Results/T_err.txt", "wb");
	FILE *err_q = fopen("Results/q_err.txt", "wb");
	FILE *err_u = fopen("Results/u_err.txt", "wb");
	FILE *err_w = fopen("Results/w_err.txt", "wb");

	// Write data to files
	for (int i = 0; i < numCellsX; ++i) {
		double x = getCellCenterX(i);
		for (int j = 0; j < numCellsP; ++j) {
			double p = getCellCenterP(i, j);

			double numer_T = T_sl[i][j], numer_q = q_sl[i][j],
					numer_u = u_sl[i][j], numer_w = w_sl[i][j];
			fprintf(res_T, "%1.20e ", numer_T);
			fprintf(res_q, "%1.20e ", numer_q);
			fprintf(res_u, "%1.20e ", numer_u);
			fprintf(res_w, "%1.20e ", numer_w);

			double exact_T = (*IC_T_fcnPtr)(x, p, finalTime),
					exact_q = (*IC_q_fcnPtr)(x, p, finalTime),
					exact_u = (*IC_u_fcnPtr)(x, p, finalTime),
					exact_w = (*IC_w_fcnPtr)(x, p, finalTime);
			fprintf(err_T, "%1.20e ", numer_T - exact_T);
			fprintf(err_q, "%1.20e ", numer_q - exact_q);
			fprintf(err_u, "%1.20e ", numer_u - exact_u);
			fprintf(err_w, "%1.20e ", numer_w - exact_w);
		}
	}

	fclose(res_T); fclose(res_q); fclose(res_u); fclose(res_w);
	fclose(err_T); fclose(err_q); fclose(err_u); fclose(err_w);
}

void printExactSolnToFile() {
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
			fprintf(exact_T, "%1.20e ", (*IC_T_fcnPtr)(x, p, t));
			fprintf(exact_q, "%1.20e ", (*IC_q_fcnPtr)(x, p, t));
			fprintf(exact_u, "%1.20e ", (*IC_u_fcnPtr)(x, p, t));
			fprintf(exact_w, "%1.20e ", (*IC_w_fcnPtr)(x, p, t));
		}
	}
	fclose(exact_T); fclose(exact_q); fclose(exact_u); fclose(exact_w);
}

// Write the result to file: res.txt for the solution and err.txt for the error
// This function is used in plotting the solution and the error
void writeResToFileForMovie_T_test(int tt) {
	if (tt % movieFrameFreq)
		return;
	char filename[20];
	sprintf(filename, "MovieFrames/T_soln_%d.txt", tt);
	FILE *res = fopen(filename, "wb");
	sprintf(filename, "MovieFrames/T_err_%d.txt", tt);
	FILE *err = fopen(filename, "wb");
	/*sprintf(filename, "MovieFrames/T_exact_%d.txt", tt);
	FILE *exact = fopen(filename, "wb");*/

	double t = tt * Dt;
	for (int i = 0; i < numCellsX; ++i) {
		double x = getCellCenterX(i);
		for (int j = 0; j < numCellsP; ++j) {
			double p = getCellCenterP(i, j);
			double numerVal = T_sl[i][j];
			fprintf(res, "%1.20e ", numerVal);

			double exactVal = (*IC_T_fcnPtr)(x, p, t);
			fprintf(err, "%1.20e ", numerVal - exactVal);
			// fprintf(exact, "%1.20e ", exactVal);

		}
	}
	fclose(res); fclose(err); //fclose(exact);
}

void printSolnToFile(int tt) {
	if (tt % movieFrameFreq)
		return;
	char filename[20];
	sprintf(filename, "MovieFrames/T_soln_%d.txt", tt);
	FILE *res_T = fopen(filename, "wb");
	sprintf(filename, "MovieFrames/q_soln_%d.txt", tt);
	FILE *res_q = fopen(filename, "wb");
	sprintf(filename, "MovieFrames/u_soln_%d.txt", tt);
	FILE *res_u = fopen(filename, "wb");
	sprintf(filename, "MovieFrames/w_soln_%d.txt", tt);
	FILE *res_w = fopen(filename, "wb");

	for (int i = 0; i < numCellsX; ++i)
		for (int j = 0; j < numCellsP; ++j) {
			fprintf(res_T, "%1.20e ", T_sl[i][j]);
			fprintf(res_q, "%1.20e ", q_sl[i][j]);
			fprintf(res_u, "%1.20e ", u_sl[i][j]);
			fprintf(res_w, "%1.20e ", w_sl[i][j]);
		}
	fclose(res_T); fclose(res_q); fclose(res_u); fclose(res_w);
}

void printExactVelocityToFile(int tt) {
	if (tt % movieFrameFreq)
		return;

	printf("\r  - Printing movie frame for tt = %d", tt);
	fflush(stdout);
	char filename[20];
	sprintf(filename, "MovieFrames/u_exact_%d.txt", tt);
	FILE *u_ex = fopen(filename, "wb");
	sprintf(filename, "MovieFrames/w_exact_%d.txt", tt);
	FILE *w_ex = fopen(filename, "wb");
	double t = tt * Dt;
	for (int i = 0; i < numCellsX; ++i) {
		double x = getCellCenterX(i);
		for (int j = 0; j < numCellsP; ++j) {
			double p = getCellCenterP(i, j);
			fprintf(u_ex, "%1.20e ", (*IC_u_fcnPtr)(x, p, t));
			fprintf(w_ex, "%1.20e ", (*IC_w_fcnPtr)(x, p, t));
		}
	}
	fclose(u_ex); fclose(w_ex);
}

void printExactVelocityToFile() {
	printf("\n>> Printing velocity field.\n");
	for (int tt = 1; tt <= numTimeStep; ++tt)
		printExactVelocityToFile(tt);
	printf("\n\nExact velocity field printed.\n");
	printParamToFile();
	printMeshToFile();
	printResToFile();
}

void peformAnalysis() {
	clock_t start = clock();

	printSchemeDescription();
	showL2Errors();
	if (printResToFile_opt) {
		printParamToFile();
		printMeshToFile();
		printResToFile();
	}
	if (printExactSolnToFile_opt)
		printExactSolnToFile();
	printf("\n- Analysis complete. Time used = %1.2fs.\n",
			((double) (clock() - start)) / CLOCKS_PER_SEC);
}

void (*aveSoln_fcnPtr)(int tt);

void aveSoln_oneTerm(double sl[numCellsX][numCellsP]) {
	for (int i = 1; i <= Nx; ++i)
		for (int j = 1; j <= Np; ++j)
			sl[i][j] = 0.5 * (sl[i][j] + sl[i - 1][j]);
}

void aveSoln(int tt) {
	if (tt % aveFreq == 0) {
		aveSoln_oneTerm(T_sl);
		// aveSoln_oneTerm(q_sl);
	}
	aveSoln_oneTerm(u_sl);
	aveSoln_oneTerm(w_sl);
	aveSoln_oneTerm(phix_sl);
}

void noAveSoln(int tt) {
	// Empty
}

void printResToFile_convAnalysis() {
	FILE *f = fopen("Results/ConvAnalysis.txt", "a");
	fprintf(f, "%d  %1.7e  %1.7e  % 1.7e  %1.7e  %1.2f\n",
			Nx, relative_L2Err_T_global, relative_L2Err_q_global,
			relative_L2Err_u_global, relative_L2Err_w_global,
			computationTime);
	fclose(f);
}

void closeFiles_analysis() {
	fclose(T_norm_filePtr); fclose(q_norm_filePtr);
	fclose(u_norm_filePtr); fclose(w_norm_filePtr);
	fclose(time_filePtr);
}

void setAnalysis() {
	T_norm_filePtr = fopen("Results/T_norm.txt", "wb");
	q_norm_filePtr = fopen("Results/q_norm.txt", "wb");
	u_norm_filePtr = fopen("Results/u_norm.txt", "wb");
	w_norm_filePtr = fopen("Results/w_norm.txt", "wb");
	time_filePtr = fopen("Results/time.txt", "wb");

	if (averageResult_opt)
		aveSoln_fcnPtr = &aveSoln;
	else
		aveSoln_fcnPtr = &noAveSoln;
}

#endif /* ANALYSIS_H_ */
