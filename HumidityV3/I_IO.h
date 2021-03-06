/*
 * Analysis.h
 *
 *  Created on: Oct 26, 2017
 *      Author: chuckjia
 *
 *  Methods and arrays to perform file and console I/O
 */

#ifndef I_IO_H_
#define I_IO_H_
#include "H_Conditions.h"

#ifdef __APPLE__
#define CURR_OS 0
#elif __linux
#define CURR_OS 1
#elif _WIN32
#define CURR_OS 2
#error "Unknown OS"
#endif

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * File I/O
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Write parameters in CSV file
void writeCSV_param() {
	FILE *f = fopen("Output/Param.csv", "wb");
	fprintf(f, "%1.20e,%1.20e,%1.20e,%1.20e,%1.20e,%d,%d,%1.20e,%d,%d,%d\n",
			x0, xf, pA, (*pB_fptr)(x0), (*pB_fptr)(xf), Nx, Np, Dt, numTimeStep, movieFrameFreq, aveSolnFreq_T);
	fclose(f);
}

void writeCSV_testParam() {
	FILE *f = fopen("Output/ParamTests.csv", "wb");
	fprintf(f, "Using compatible IC for q,");
	if (_useCompatibleInitQ_) fprintf(f, "true"); else fprintf(f, "false");
	fprintf(f, "\n");
	fprintf(f, "Enforcing top BC for u and w,");
	if (_enforceTopBC_) fprintf(f, "true"); else fprintf(f, "false");
	fprintf(f, "\n");

	fprintf(f, "aveSolnFreq_T,%d\n",    aveSolnFreq_T);
	fprintf(f, "aveSolnFreq_q,%d\n",    aveSolnFreq_q);
	fprintf(f, "aveSolnFreq_u,%d\n",    aveSolnFreq_u);
	fprintf(f, "aveSolnFreq_w,%d\n",    aveSolnFreq_w);
	fprintf(f, "aveSolnFreq_phix,%d\n", aveSolnFreq_phix);

	fprintf(f, "numMountains,%d\n", numMountain);
	fprintf(f, "Mountain 1 Height,%1.20e\n", _h_pB_MDL0);
	fprintf(f, "Mountain 2 Height,%1.20e\n", _h2_pB_MDL0);
	fclose(f);
}

// Print the mesh grid points to CSV file
// The x- and p-coordinates of the grid points are stored in two different files, with names filenameX and filenameP.
void writeCSV_meshGrid(const char* filenameX, const char* filenameP) {
	FILE *fx = fopen(filenameX, "wb");
	FILE *fp = fopen(filenameP, "wb");

	int last_j = numGridPtP - 1;  // Used to print the last entry of each row, which does not have commas following

	for (int i = 0; i < numGridPtX; ++i) {
		double gridX = meshGridX_[i];
		for (int j = 0; j < last_j; ++j) {
			double gridP = meshGridP_[i][j];
			fprintf(fx, "%1.20e,", gridX);
			fprintf(fp, "%1.20e,", gridP);
		}
		fprintf(fx, "%1.20e\n", gridX);
		fprintf(fp, "%1.20e\n", meshGridP_[i][last_j]);
	}

	fclose(fx); fclose(fp);
}

// Print the mesh grid points to files with standard names
void writeCSV_meshGrid() {
	writeCSV_meshGrid("Output/MeshGrid_X.csv", "Output/MeshGrid_P.csv");
}

// Print cell center coordinates to file
void writeCSV_cellCenters() {
	// Print x-coordinates
	FILE *f = fopen("Output/CellCenters_X.csv", "wb");
	int last_j = numCellP - 1;
	for (int i = 0; i < numCellX; ++i) {
		double centerX = cellCenterX_[i];
		for (int j = 0; j < last_j; ++j)
			fprintf(f, "%1.20e,", centerX);
		fprintf(f, "%1.20e\n", centerX);
	}
	fclose(f);

	// Print p-coordinates
	writeCSV_matrix(cellCenterP_, "Output/CellCenters_P.csv");
}

// Print solutions to movie frame files. File names are appended with time step numbers
void writeMovie_soln(int tt) {
	string suffix = "_" + to_string(tt) + ".csv";
	string filename;

	filename = movieFramesFolderName + "/T" + suffix;
	writeCSV_matrix(T_, filename.c_str());

	filename = movieFramesFolderName + "/q" + suffix;
	writeCSV_matrix(q_, filename.c_str());

	filename = movieFramesFolderName + "/u" + suffix;
	writeCSV_matrix(u_, filename.c_str());

	filename = movieFramesFolderName + "/w" + suffix;
	writeCSV_matrix(w_, filename.c_str());
}


/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Console I/O
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Print title messages
void printTitle() {
	printf("* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== \n");
	printf("* Solving Atmospherical Model Using Upwind-type Godunov Scheme\n");
	printf("* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== \n\n");
}

// Print all parameters for diagnostics
void printDiagnostics() {
	printf("\n>> DIAGNOSTIC INFORMATION\n\nParameters: \n");

	printf("\n - Time\n");
	printf("    [1] Dt = %1.2e, finalTime = %1.2es, numTimeStep = %d\n", Dt, finalTime, numTimeStep);

	printf("\n - Domain geometry\n");
	printf("    [1] x0 = %1.2f, xf = %1.2f \n", x0, xf);
	printf("    [2] pA = %1.2f, pB(x0) = %1.2f, pB[(x0+xf)/2] = %1.2f, pB(xf) = %1.2f\n",
			pA, (*pB_fptr)(x0), (*pB_fptr)(0.5 * (x0 + xf)), (*pB_fptr)(xf));

	printf("\n - Mesh specifications\n");
	printf("    [1] Nx = %d, Np = %d\n", Nx, Np);
	printf("    [2] Dx = %1.2f\n", Dx);
	printf("    [3] Dp(x0) = %1.2f, Dp[(x0+xf)/2] = %1.2f, Dp(xf) = %1.2f\n",
			getCellLeftDp(1), getCellLeftDp(lastRealIndexX / 2), getCellLeftDp(lastRealIndexX));
}

// Print a summary of the scheme description
void printSchemeSummary() {
	printf("\n>> SCHEME DESCRIPTION\n");
	printf("\n  MODEL %d:  \n", modelNo);

	printf("  [1] Time           : ");
	printf("Dt = %1.2e,  finalTime = %1.2es,  numTimeSteps = %d\n", Dt, finalTime, numTimeStep);

	printf("  [2] Mesh Specs     : ");
	printf("Nx = %3d,       Np = %3d\n", Nx, Np);

	printf("  [3] Average Method : ");
	if (aveMethodApplied) {
		printf("every %3d, %3d, %3d, %3d, %3d steps on T, q, u, w, phi_x, resp.\n", aveSolnFreq_T, aveSolnFreq_q, aveSolnFreq_u, aveSolnFreq_w, aveSolnFreq_phix);
		printf("                               |    |    |    |    |\n");
		printf("                               T    q    u    w  phi_x\n");
	} else
		printf("No averaging is applied.\n");
}


/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Calculate and Print Errors/Exact Solutions: For Test Models Only
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Calculate the numerical L2 error with in one cell by comparing with the exact solution
double absL2ErrInACell_helper(double numerVal, double exactVal, double vol) {
	double err = exactVal - numerVal;
	return vol * err * err;
}

const double GRTD_PREC = 1e-15;  // Guaranteed precision

// Helper for computing the relative L2 error. This helper ensures if the numerator and the denominator are both too small,
// then an relative error of 0 is returned
double relatL2Err_helper(double num, double denom) {
	if (denom > GRTD_PREC)
		return sqrt(num / denom);
	else if (num < GRTD_PREC)
		return 0;
	return -987654321e10;
}

double relatL2Err_T_, relatL2Err_q_, relatL2Err_u_, relatL2Err_w_;

// Print the L2 relative and absolute errors
// !!AlphaVersion!!
void showL2Errors(double t) {
	double num_T = 0, num_q = 0, num_u = 0, num_w = 0,
			denom_T = 0, denom_q = 0, denom_u = 0, denom_w = 0;
	for (int i = 1; i < Nx; ++i) {
		double x = getCellCenterX(i), vol = getCellVol(i);
		for (int j = 1; j < Np; ++j) {
			double p = getCellCenterP(i, j);
			double exactVal_T = (*initT_fptr)(x, p, t),
					exactVal_q = (*initQ_fptr)(x, p, t),
					exactVal_u = (*initU_fptr)(x, p, t),
					exactVal_w = (*initW_fptr)(x, p, t);
			// Calculate errors
			num_T += absL2ErrInACell_helper(exactVal_T, T_[i][j], vol);
			num_q += absL2ErrInACell_helper(exactVal_q, q_[i][j], vol);
			num_u += absL2ErrInACell_helper(exactVal_u, u_[i][j], vol);
			num_w += absL2ErrInACell_helper(exactVal_w, w_[i][j], vol);
			denom_T += vol * exactVal_T * exactVal_T;
			denom_q += vol * exactVal_q * exactVal_q;
			denom_u += vol * exactVal_u * exactVal_u;
			denom_w += vol * exactVal_w * exactVal_w;
		}
	}
	double absL2Err_T = sqrt(num_T), absL2Err_q = sqrt(num_q),
			absL2Err_u = sqrt(num_u), absL2Err_w = sqrt(num_w);
	double relatL2Err_T = relatL2Err_helper(num_T, denom_T),
			relatL2Err_q = relatL2Err_helper(num_q, denom_q),
			relatL2Err_u = relatL2Err_helper(num_u, denom_u),
			relatL2Err_w = relatL2Err_helper(num_w, denom_w);

	printf("\n  - L2 relative error for T, q, u, w = [%1.7e, %1.7e, %1.7e, %1.7e]\n",
			relatL2Err_T, relatL2Err_q, relatL2Err_u, relatL2Err_w);
	printf("  - L2 discrete norm for T, q, u, w = [%1.7e, %1.7e, %1.7e, %1.7e]\n",
			sqrt(denom_T), sqrt(denom_q), sqrt(denom_u), sqrt(denom_w));
	printf("  - L2 absolute error for T, q, u, w = [%1.7e, %1.7e, %1.7e, %1.7e]\n",
			absL2Err_T, absL2Err_q, absL2Err_u, absL2Err_w);
	relatL2Err_T_ = relatL2Err_T;
	relatL2Err_q_ = relatL2Err_q;
	relatL2Err_u_ = relatL2Err_u;
	relatL2Err_w_ = relatL2Err_w;
}

// Print to console L2 errors at the final time of the computation
void showL2Errors() {
	printf("\n>> NUMERICAL ERRORS\n");
	showL2Errors(finalTime);
}

void writeCSV_realTimeErr(int tt) {
	string suffix = "_" + to_string(tt) + ".csv";
	string filename;
	double t = Dt * tt;

	// Numerical errors files
	filename = movieFramesFolderName + "/err_T" + suffix;
	FILE *err_T = fopen(filename.c_str(), "wb");
	filename = movieFramesFolderName + "/err_q" + suffix;
	FILE *err_q = fopen(filename.c_str(), "wb");
	filename = movieFramesFolderName + "/err_u" + suffix;
	FILE *err_u = fopen(filename.c_str(), "wb");
	filename = movieFramesFolderName + "/err_w" + suffix;
	FILE *err_w = fopen(filename.c_str(), "wb");

	// Write data to files
	int last_j = numCellP - 1;
	for (int i = 0; i < numCellX; ++i) {
		double x = getCellCenterX(i);
		for (int j = 0; j < last_j; ++j) {
			double p = getCellCenterP(i, j),
					numer_T = T_[i][j],
					numer_q = q_[i][j],
					numer_u = u_[i][j],
					numer_w = w_[i][j],
					exact_T = (*initT_fptr)(x, p, t),
					exact_q = (*initQ_fptr)(x, p, t),
					exact_u = (*initU_fptr)(x, p, t),
					exact_w = (*initW_fptr)(x, p, t);
			fprintf(err_T, "%1.20e,", numer_T - exact_T);
			fprintf(err_q, "%1.20e,", numer_q - exact_q);
			fprintf(err_u, "%1.20e,", numer_u - exact_u);
			fprintf(err_w, "%1.20e,", numer_w - exact_w);
		}  // End of loop for j

		double p = getCellCenterP(i, last_j),
				numer_T = T_[i][last_j],
				numer_q = q_[i][last_j],
				numer_u = u_[i][last_j],
				numer_w = w_[i][last_j],
				exact_T = (*initT_fptr)(x, p, t),
				exact_q = (*initQ_fptr)(x, p, t),
				exact_u = (*initU_fptr)(x, p, t),
				exact_w = (*initW_fptr)(x, p, t);
		fprintf(err_T, "%1.20e\n", numer_T - exact_T);
		fprintf(err_q, "%1.20e\n", numer_q - exact_q);
		fprintf(err_u, "%1.20e\n", numer_u - exact_u);
		fprintf(err_w, "%1.20e\n", numer_w - exact_w);
	}
	printf("\n>> Writing real time error as CSV files.\n");
	fclose(err_T); fclose(err_q); fclose(err_u); fclose(err_w);
}

// Write the result to file: res.txt for the solution and err.txt for the error
// This function is used in plotting the solution and the error
// !!AlphaVersion!!
void writeCSV_finalSolnErr() {
	// Numerical solution files
	FILE *res_T = fopen("Output/T_soln.csv", "wb");
	FILE *res_q = fopen("Output/q_soln.csv", "wb");
	FILE *res_u = fopen("Output/u_soln.csv", "wb");
	FILE *res_w = fopen("Output/w_soln.csv", "wb");
	// Numerical errors files
	FILE *err_T = fopen("Output/T_err.csv", "wb");
	FILE *err_q = fopen("Output/q_err.csv", "wb");
	FILE *err_u = fopen("Output/u_err.csv", "wb");
	FILE *err_w = fopen("Output/w_err.csv", "wb");

	// Write data to files
	int last_j = numCellP - 1;
	for (int i = 0; i < numCellX; ++i) {
		double x = getCellCenterX(i);
		for (int j = 0; j < last_j; ++j) {
			double p = getCellCenterP(i, j),
					numer_T = T_[i][j], numer_q = q_[i][j], numer_u = u_[i][j], numer_w = w_[i][j];
			fprintf(res_T, "%1.20e,", numer_T);
			fprintf(res_q, "%1.20e,", numer_q);
			fprintf(res_u, "%1.20e,", numer_u);
			fprintf(res_w, "%1.20e,", numer_w);

			double exact_T = (*initT_fptr)(x, p, finalTime),
					exact_q = (*initQ_fptr)(x, p, finalTime),
					exact_u = (*initU_fptr)(x, p, finalTime),
					exact_w = (*initW_fptr)(x, p, finalTime);
			fprintf(err_T, "%1.20e,", numer_T - exact_T);
			fprintf(err_q, "%1.20e,", numer_q - exact_q);
			fprintf(err_u, "%1.20e,", numer_u - exact_u);
			fprintf(err_w, "%1.20e,", numer_w - exact_w);
		}

		double p = getCellCenterP(i, last_j),
				numer_T = T_[i][last_j], numer_q = q_[i][last_j], numer_u = u_[i][last_j], numer_w = w_[i][last_j];
		fprintf(res_T, "%1.20e\n", numer_T);
		fprintf(res_q, "%1.20e\n", numer_q);
		fprintf(res_u, "%1.20e\n", numer_u);
		fprintf(res_w, "%1.20e\n", numer_w);

		double exact_T = (*initT_fptr)(x, p, finalTime),
				exact_q = (*initQ_fptr)(x, p, finalTime),
				exact_u = (*initU_fptr)(x, p, finalTime),
				exact_w = (*initW_fptr)(x, p, finalTime);
		fprintf(err_T, "%1.20e\n", numer_T - exact_T);
		fprintf(err_q, "%1.20e\n", numer_q - exact_q);
		fprintf(err_u, "%1.20e\n", numer_u - exact_u);
		fprintf(err_w, "%1.20e\n", numer_w - exact_w);
	}

	fclose(res_T); fclose(res_q); fclose(res_u); fclose(res_w);
	fclose(err_T); fclose(err_q); fclose(err_u); fclose(err_w);
}

// Print the exact solution to file at final time
// !!AlphaVersion!!
void writeCSV_exactSoln() {
	// Write final numerical solution to files
	FILE *exact_T = fopen("Output/T_exact.csv", "wb");
	FILE *exact_q = fopen("Output/q_exact.csv", "wb");
	FILE *exact_u = fopen("Output/u_exact.csv", "wb");
	FILE *exact_w = fopen("Output/w_exact.csv", "wb");

	int last_j = numCellP - 1;
	for (int i = 0; i < numCellX; ++i) {
		double x = getCellCenterX(i);
		for (int j = 0; j < last_j; ++j) {
			double p = getCellCenterP(i, j);
			fprintf(exact_T, "%1.20e,", (*initT_fptr)(x, p, finalTime));
			fprintf(exact_q, "%1.20e,", (*initQ_fptr)(x, p, finalTime));
			fprintf(exact_u, "%1.20e,", (*initU_fptr)(x, p, finalTime));
			fprintf(exact_w, "%1.20e,", (*initW_fptr)(x, p, finalTime));
		}

		double p = getCellCenterP(i, last_j);
		fprintf(exact_T, "%1.20e\n", (*initT_fptr)(x, p, finalTime));
		fprintf(exact_q, "%1.20e\n", (*initQ_fptr)(x, p, finalTime));
		fprintf(exact_u, "%1.20e\n", (*initU_fptr)(x, p, finalTime));
		fprintf(exact_w, "%1.20e\n", (*initW_fptr)(x, p, finalTime));
	}

	fclose(exact_T); fclose(exact_q); fclose(exact_u); fclose(exact_w);
}


/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Calculate and Print Norms of Solutions: For Physical Models
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Print the L2 norm of a solution
double calcL2Norm(double sl[numCellX][numCellP]) {
	double norm = 0;
	for (int i = 1; i <= Nx; ++i) {
		double vol = getCellVol(i);
		for (int j = 1; j <= Np; ++j) {
			double slVal = sl[i][j];
			norm += vol * slVal * slVal;
		}
	}
	return sqrt(norm);
}

FILE *T_norm_file, *q_norm_file, *u_norm_file, *w_norm_file, *time_file;

// Calculate L2 norms of all solutions and print to file
// !!AlphaVersion!!
void calcL2Norm(double t) {
	double norm_T = 0, norm_q = 0, norm_u = 0, norm_w = 0;
	for (int i = 1; i <= Nx; ++i) {
		double vol = getCellVol(i);
		for (int j = 1; j <= Np; ++j) {
			double slVal;
			slVal = T_[i][j];
			norm_T += vol * slVal * slVal;
			slVal = q_[i][j];
			norm_q += vol * slVal * slVal;
			slVal = u_[i][j];
			norm_u += vol * slVal * slVal;
			slVal = w_[i][j];
			norm_w += vol * slVal * slVal;
		}
	}
	norm_T = sqrt(norm_T);
	norm_q = sqrt(norm_q);
	norm_u = sqrt(norm_u);
	norm_w = sqrt(norm_w);

	printf("      - L2 norm (T, q, u, w) = (%1.2f, %1.2f, %1.2f, %1.2f)", norm_T - 2.055e6, norm_q, norm_u, norm_w);
	fprintf(T_norm_file, "%1.20e\n", norm_T);
	fprintf(q_norm_file, "%1.20e\n", norm_q);
	fprintf(u_norm_file, "%1.20e\n", norm_u);
	fprintf(w_norm_file, "%1.20e\n", norm_w);
	fprintf(time_file, "%1.20e\n", t);
}

double computationTime = 0;

// !!AlphaVersion!!
void printExactVelocityToMovie(int tt) {
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
	for (int i = 0; i < numCellX; ++i) {
		double x = getCellCenterX(i);
		for (int j = 0; j < numCellP; ++j) {
			double p = getCellCenterP(i, j);
			fprintf(u_ex, "%1.20e ", (*initU_fptr)(x, p, t));
			fprintf(w_ex, "%1.20e ", (*initW_fptr)(x, p, t));
		}
	}
	fclose(u_ex); fclose(w_ex);
}

// !!AlphaVersion!!
void printExactVelocityToFile() {
	printf("\n>> Printing velocity field.\n");
	for (int tt = 1; tt <= numTimeStep; ++tt)
		printExactVelocityToMovie(tt);
	printf("\n\nExact velocity field printed.\n");
	writeCSV_param();
	writeCSV_finalSolnErr();
}

// !!AlphaVersion!!
void peformAnalysis() {
	clock_t start = clock();

	printSchemeSummary();
	showL2Errors();
	if (_printResultToFile_) {
		writeCSV_param();
		writeCSV_finalSolnErr();
	}
	if (_printExactSolnToFile_)
		writeCSV_exactSoln();
	printf("\n- Analysis complete. Time used = %1.2fs.\n", ((double) (clock() - start)) / CLOCKS_PER_SEC);
}

// !!AlphaVersion!!
void writeCSV_convAnalysis() {
	FILE *f = fopen("Results/ConvAnalysis.csv", "a");
	fprintf(f, "%d,%1.7e,%1.7e,%1.7e,%1.7e,%1.7e\n",
			Nx, relatL2Err_T_, relatL2Err_q_,
			relatL2Err_u_, relatL2Err_w_,
			computationTime);
	fclose(f);
}

// !!AlphaVersion!!
void closeGlobalFiles_IO() {
	fclose(T_norm_file); fclose(q_norm_file);
	fclose(u_norm_file); fclose(w_norm_file);
	fclose(time_file);
}

void printTestValue(double mat[Nx + 2][Np + 2], int i, int j, int target_tt, int tt) {
	if (target_tt != tt)
		return;
	printf("\nCurrent time iteration = %d\n", tt);
	if (mat == T_)         printf("T");
	else if (mat == q_)    printf("q");
	else if (mat == u_)    printf("u");
	else if (mat == w_)    printf("w");
	else if (mat == phix_) printf("phix");
	printf("[%d][%d] = %1.25e\n", i, j, mat[i][j]);
}

void printTestValue(double mat[Nx + 2][Np + 2], int i, int j, int target_tt, int tt, string msg) {
	printf("%s", msg.c_str());
	printTestValue(mat, i, j, target_tt, tt);
}

void printTestValue(double mat[Nx + 2][Np + 2], int i, int j) {
	printTestValue(mat, i, j, -1, -1);
}

void printTestValue(double mat[Nx + 2][Np + 2], int i, int j, string msg) {
	printf("%s", msg.c_str());
	printTestValue(mat, i, j, -1, -1);
}

// !!AlphaVersion!!
void setAnalysis() {
	T_norm_file = fopen("Output/T_norm.csv", "wb");
	q_norm_file = fopen("Output/q_norm.csv", "wb");
	u_norm_file = fopen("Output/u_norm.csv", "wb");
	w_norm_file = fopen("Output/w_norm.csv", "wb");
	time_file   = fopen("Output/time.csv",   "wb");
}

#endif /* I_IO_H_ */
