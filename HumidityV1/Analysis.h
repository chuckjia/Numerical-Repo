/*
 * Analysis.h
 *
 *  Created on: Oct 26, 2017
 *      Author: chuckjia
 */

#ifndef ANALYSIS_H_
#define ANALYSIS_H_
#include "TimeSteps.h"

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Print all parameters for diagnostics
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void printDiagnostics() {
	printf(">> DIAGNOSTIC INFORMATION\n");
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

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Calculate and Print Errors
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

double cellErr_helper(int i, int j, double exactVal, double vol) {
	return vol * pow(exactVal - T_sl[i][j], 2);
}

double relativeL2Err_helper(double num, double denom) {
	if (denom > 1e-15)
		return sqrt(num / denom);
	else if (num < 1e-15)
		return 0;
	return -99;
}

// Print the L2 relative and absolute errors
void showL2Errors() {
	double t = finalTime;
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
			num_T += cellErr_helper(i, j, exactVal_T, vol);
			num_q += cellErr_helper(i, j, exactVal_q, vol);
			num_u += cellErr_helper(i, j, exactVal_u, vol);
			num_w += cellErr_helper(i, j, exactVal_w, vol);
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
	printf("\n- L2 relative error for T, q, u = [%1.4e, %1.4e, %1.4e, %1.4e]\n",
			relativeL2Err_T, relativeL2Err_q, relativeL2Err_u, relativeL2Err_w);
	printf("\n- L2 absolute error for T, q, u = [%1.4e, %1.4e, %1.4e, %1.4e]\n",
			absL2Err_T, absL2Err_q, absL2Err_u, absL2Err_w);
}

void printParToFile() {
	FILE *f = fopen("Results/par.txt", "wb");
	fprintf(f, "%1.15f %1.15f %1.15f %1.15f %1.15f %d %d %1.15f %d",
			x0, xf, pA, (*pB_fcnPtr)(x0), (*pB_fcnPtr)(xf), Nx, Np, Dt, numTimeSteps);
	fclose(f);
}

void writeCellCentersToFile() {
	FILE *fCellCenterX = fopen("Results/cellCentersX.txt", "wb");
	FILE *fCellCenterP = fopen("Results/cellCentersP.txt", "wb");
	for (int i = 1; i <= Nx; ++i)
		for (int j = 1; j <= Np; ++j) {
			fprintf(fCellCenterX, "%f ", getCellCenterX(i, j));
			fprintf(fCellCenterP, "%f ", getCellCenterP(i, j));
		}
	fclose(fCellCenterX); fclose(fCellCenterP);
}

// Write the result to file: res.txt for the solution and err.txt for the error
// This function is used in ploting the solution and the error
void writeResToFile() {
	// Write final numerical solution to files
	FILE *res_T = fopen("Results/res_T.txt", "wb");
	FILE *res_q = fopen("Results/res_q.txt", "wb");
	FILE *res_u = fopen("Results/res_u.txt", "wb");
	FILE *res_w = fopen("Results/res_w.txt", "wb");
	for (int i = 1; i < lastGhostIndexX; ++i)
		for (int j = 1; j < lastGhostIndexP; ++j) {
			fprintf(res_T, "%f ", T_sl[i][j]);
			fprintf(res_q, "%f ", q_sl[i][j]);
			fprintf(res_u, "%f ", u_sl[i][j]);
			fprintf(res_w, "%f ", w_sl[i][j]);
		}
	fclose(res_T); fclose(res_q); fclose(res_u); fclose(res_w);

	// Write final numerical errors to files
	FILE *err_T = fopen("Results/err_T.txt", "wb");
	FILE *err_q = fopen("Results/err_q.txt", "wb");
	FILE *err_u = fopen("Results/err_u.txt", "wb");
	FILE *err_w = fopen("Results/err_w.txt", "wb");
	for (int i = 1; i < lastGhostIndexX; ++i) {
		double x = getCellCenterX(i);
		for (int j = 1; j < lastGhostIndexP; ++j) {
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
	printDiagnostics();
	showL2Errors();
	printParToFile();
	writeCellCentersToFile();
	writeResToFile();
}

#endif /* ANALYSIS_H_ */
