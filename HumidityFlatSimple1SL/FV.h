/*
 * FV.h
 *
 *  Created on: Jul 18, 2017
 *      Author: chuckjia
 */

#ifndef FV_H_
#define FV_H_

#include <stdio.h>
#include <math.h>

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * General Settings
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

/* ----- ----- ----- ----- ----- -----
 * Scheme specifications
 * ----- ----- ----- ----- ----- ----- */

// Number of space steps
const int numDivisions = 200;
// Number of time steps
const int numTimeSteps = 200;
// Size of time steps
const double Dt = 1e-5;
const double finalTime = numTimeSteps * Dt;

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Other Constants
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

const double x0 = 0;
const double xL = 10;
const double pA = 0;
const double pB = 10;

const int Nx = numDivisions;
const int Np = numDivisions;
const int numCellsX = Nx + 2, numCellsP = Np + 2;
const int lastRealIndexX = Nx, lastRealIndexP = Np;
const int lastGhostIndexX = Nx + 1, lastGhostIndexP = Np + 1;

const double Dx = (xL - x0) / Nx, Dp = (pB - pA) / Np;
const double cellVol = Dx * Dp;

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Getter Functions For the Mesh
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Return the x coordinate of cell left side
double getCellLeftX(int i, int j) {
	return x0 + (i - 1) * Dx;
}

// Return the x coordinate of cell right side
double getCellRightX(int i, int j) {
	return x0 + i * Dx;
}

// Return the p coordinate of cell bottom side
double getCellBottP(int i, int j) {
	return pA + (j - 1) * Dp;
}

// Return the p coordinate of cell top side
double getCellTopP(int i, int j) {
	return pA + j * Dp;
}

// Return the x coordinate of cell center
double getCellCenterX(int i, int j) {
	return x0 + (i - 0.5) * Dx;
}

// Return the p coordinate of cell center
double getCellCenterP(int i, int j) {
	return pA + (j - 0.5) * Dp;
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Model Functions
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

const double c1_uFcn = M_PI / pB, c2_uFcn = 2 * M_PI / xL;

double u_fcn(double x, double p) {
	return 7.5 + cos(c1_uFcn * p) * cos(c2_uFcn * x);
}

double u_xDer_fcn(double x, double p) {
	return -cos(c1_uFcn * p) * c2_uFcn * sin(c2_uFcn * x);
}

const double c1_wFcn = 2 * M_PI / pB, c2_wFcn = 2 * M_PI / xL;

double w_fcn(double x, double p) {
	return sin(c1_wFcn * p) * cos(c2_wFcn * x);
}

double w_pDer_fcn(double x, double p) {
	return c1_wFcn * cos(c1_wFcn * p) * cos(c2_wFcn * x);
}

const double c1_Test1 = 0.2 * 2 * M_PI;
const double c2_Test1 = 0.2 * 2 * M_PI;
const double c3_Test1 = 2 * M_PI;

double exact_T_Test1(double x, double p, double t) {
	return sin(c1_Test1 * x) * sin(c2_Test1 * p) * cos(c3_Test1 * t);
}

double boundaryVal = 0;

double source_Test1(double T, double x, double p, double t) {
	double xInput = c1_Test1 * x, pInput = c2_Test1 * p, tInput = c3_Test1 * t;
	double xPart = sin(xInput), pPart = sin(pInput), tPart = cos(tInput);
	return - xPart * pPart * c3_Test1 * sin(tInput) +
			T * (u_xDer_fcn(x, p) + w_pDer_fcn(x, p)) +
			c1_Test1 * cos(c1_Test1 * x) * pPart * tPart * u_fcn(x, p) +
			xPart * c2_Test1 * cos(c2_Test1 * p) * tPart * w_fcn(x, p);
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Initial Conditions
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

double sl[numCellsX][numCellsP];

void enforceInitCond() {
	for (int i = 0; i < numCellsX; i++)
		for (int j = 0; j < numCellsP; j++) {
			double x = getCellCenterX(i, j), p = getCellCenterP(i, j);
			sl[i][j] = exact_T_Test1(x, p, 0);
		}
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Fluxes
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

double FF_Left[numCellsX][numCellsP], FF_Right[numCellsX][numCellsP];
double GG_Bott[numCellsX][numCellsP], GG_Top[numCellsX][numCellsP];

void calcFluxes_Dirichlet() {
	for (int i = 1; i < lastGhostIndexX; i++)
		for (int j = 1; j < lastGhostIndexP; j++) {
			double x, p;
			// GG fluxes: bottom
			x = getCellCenterX(i, j);
			p = getCellBottP(i, j);
			if (j == 1)
				GG_Bott[i][j] = -Dx * w_fcn(x, p) * boundaryVal;
			else
				GG_Bott[i][j] = -Dx * w_fcn(x, p) * 0.5 * (sl[i][j - 1] + sl[i][j]);

			// GG fluxes: top
			p = getCellTopP(i, j);
			if (j == lastRealIndexP)
				GG_Top[i][j] = Dx * w_fcn(x, p) * boundaryVal;
			else
				GG_Top[i][j] = Dx * w_fcn(x, p) * 0.5 * (sl[i][j] + sl[i][j + 1]);

			// FF fluxes: left
			p = getCellCenterP(i, j);
			x = getCellLeftX(i, j);
			if (i == 1)
				FF_Left[i][j] = -Dp * u_fcn(x, p) * boundaryVal;
			else
				FF_Left[i][j] = -Dp * u_fcn(x, p) * 0.5 * (sl[i - 1][j] + sl[i][j]);

			// FF fluxes: right
			x = getCellRightX(i, j);
			if (i == lastRealIndexX)
				FF_Right[i][j] = Dp * u_fcn(x, p) * boundaryVal;
			else
				FF_Right[i][j] = Dp * u_fcn(x, p) * 0.5 * (sl[i][j] + sl[i + 1][j]);
		}
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Time Steps
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void forwardEuler() {
	enforceInitCond();
	for (int tt = 0; tt < numTimeSteps; tt++) {
		double t = tt * Dt;
		//printf("\rProgress: %1.1f%%", t / finalTime * 100);
		calcFluxes_Dirichlet();
		for (int i = 1; i < lastGhostIndexX; i++)
			for (int j = 1; j < lastGhostIndexP; j++) {
				double T = sl[i][j], x = getCellCenterX(i, j), p = getCellCenterP(i, j);
				double RHS = source_Test1(T, x, p, t) -
						(GG_Bott[i][j] + GG_Top[i][j] + FF_Left[i][j] + FF_Right[i][j]) / cellVol;
				sl[i][j] += Dt * RHS;
			}
	}
	printf("\nCompleted: 100%%\n");
}

void showL2Errors() {
	double t = finalTime;
	double sum1 = 0, sum2 = 0;
	for (int i = 1; i < lastGhostIndexX; i++)
		for (int j = 1; j < lastGhostIndexP; j++) {
			double x = getCellCenterX(i, j), p = getCellCenterP(i, j);
			double ExactVal = exact_T_Test1(x, p, t);
			double NumericalVal = sl[i][j];
			sum1 += pow(ExactVal - NumericalVal, 2);
			sum2 += ExactVal * ExactVal;
		}
	double absError = sqrt(cellVol * sum1);
	double relativeError = 0;
	if (sum1 > 1e-15)
		relativeError = sqrt(sum1 / sum2);
	printf("\n- L2 relative error = %1.15f\n", relativeError);
	printf("\n- L2 absolute error = %1.15f\n", absError);
}

void printMsg() {
	printf("\nnumDivisions = %d, numTimeSteps = %d, Dt = %f, finalTime = %1.4f\n",
			numDivisions, numTimeSteps, Dt, finalTime);
}

void writeResToFile() {
	FILE *f = fopen("res.txt", "wb");
	for (int i = 1; i < lastGhostIndexX; i++)
		for (int j = 1; j < lastGhostIndexP; j++)
			fprintf(f, "%f ", sl[i][j]);
	fclose(f);
	FILE *g = fopen("err.txt", "wb");
	for (int i = 1; i < lastGhostIndexX; i++)
		for (int j = 1; j < lastGhostIndexP; j++) {
			double x = getCellCenterX(i, j), p = getCellCenterP(i, j);
			double exactSoln = exact_T_Test1(x, p, finalTime);
			fprintf(g, "%f ", sl[i][j] - exactSoln);
		}
	fclose(g);
}

#endif /* FV_H_ */
