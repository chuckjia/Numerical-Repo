/*
 * FV.h
 *
 *  Created on: Jul 18, 2017
 *      Author: chuckjia
 */

#ifndef FV_H_
#define FV_H_
#include "Mesh.h"

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

double FF[numCellsX][numCellsP];
double GG[numCellsX][numCellsP];

void calcFluxes_ClassFV_Dirichlet() {
	/*
	 * GG fluxes
	 */
	for (int i = 1; i < lastGhostIndexX; i++)
		for (int j = 1; j < lastRealIndexP; j++) {
			double x = getCellCenterX(i, j), p = getCellTopP(i, j);
			GG[i][j] = Dx * w_fcn(x, p) * 0.5 * (sl[i][j] + sl[i][j + 1]);
		}
	// Top side and bottom side of domain
	int jRange[2] = {0, lastRealIndexP};
	for (int jj = 0; jj < 2; jj++) {
		int j = jRange[jj];
		for (int i = 1; i < lastGhostIndexX; i++) {
			double x = getCellCenterX(i, j), p = getCellTopP(i, j);
			GG[i][j] = Dx * w_fcn(x, p) * boundaryVal;
		}
	}

	/*
	 * FF fluxes
	 */
	for (int i = 1; i < lastRealIndexX; i++)
		for (int j = 1; j < lastGhostIndexP; j++) {
			double x = getCellRightX(i, j), p = getCellCenterP(i, j);
			FF[i][j] = Dp * u_fcn(x, p) * 0.5 * (sl[i][j] + sl[i + 1][j]);
		}
	// Left and right sides of domain
	int iRange[2] = {0, lastRealIndexX};
	for (int ii = 0; ii < 2; ii++) {
		int i = iRange[ii];
		for (int j = 1; j < lastGhostIndexP; j++) {
			double x = getCellRightX(i, j), p = getCellCenterP(i, j);
			FF[i][j] = Dp * u_fcn(x, p) * boundaryVal;
		}
	}
}

void calcFluxes_Godunov_Dirichlet() {
	/*
	 * GG fluxes
	 */
	for (int i = 1; i < lastGhostIndexX; i++)
		for (int j = 1; j < lastRealIndexP; j++) {
			double x = getCellCenterX(i, j), p = getCellTopP(i, j);
			double wVal = w_fcn(x, p);
			if (wVal >= 0)
				GG[i][j] = Dx * wVal * sl[i][j];
			else
				GG[i][j] = Dx * wVal * sl[i][j + 1];
		}
	// Top side and bottom side of domain
	int jRange[2] = {0, lastRealIndexP};
	for (int jj = 0; jj < 2; jj++) {
		int j = jRange[jj];
		for (int i = 1; i < lastGhostIndexX; i++) {
			double x = getCellCenterX(i, j), p = getCellTopP(i, j);
			GG[i][j] = Dx * w_fcn(x, p) * boundaryVal;
		}
	}

	/*
	 * FF fluxes
	 */
	for (int i = 1; i < lastRealIndexX; i++)
		for (int j = 1; j < lastGhostIndexP; j++) {
			double x = getCellRightX(i, j), p = getCellCenterP(i, j);
			double uVal = u_fcn(x, p);
			if (uVal >= 0)
				FF[i][j] = Dp * uVal * sl[i][j];
			else
				FF[i][j] = Dp * uVal * sl[i + 1][j];
		}
	// Left and right sides of domain
	int iRange[2] = {0, lastRealIndexX};
	for (int ii = 0; ii < 2; ii++) {
		int i = iRange[ii];
		for (int j = 1; j < lastGhostIndexP; j++) {
			double x = getCellRightX(i, j), p = getCellCenterP(i, j);
			FF[i][j] = Dp * u_fcn(x, p) * boundaryVal;
		}
	}
}

double GG_Bott[numCellsX][numCellsP], GG_Top[numCellsX][numCellsP];
double FF_Left[numCellsX][numCellsP], FF_Right[numCellsX][numCellsP];

void calcFluxes_Godunov_Dirichlet2() {
	/*
	 * GG fluxes
	 */
	for (int i = 1; i < lastGhostIndexX; i++)
		for (int j = 1; j < lastGhostIndexP; j++) {
			double x, p, uVal, wVal;

			/*
			 * GG fluxes
			 */
			x = getCellCenterX(i, j);

			// GG bottom
			p = getCellBottP(i, j);
			wVal = w_fcn(x, p);
			if (j == 1)
				GG_Bott[i][j] = -Dx * wVal * boundaryVal;
			else if (wVal >= 0) {
				GG_Bott[i][j] = -Dx * wVal * sl[i][j - 1];
			}
			else {
				GG_Bott[i][j] = -Dx * wVal * sl[i][j];
			}
			// GG top
			p = getCellTopP(i, j);
			wVal = w_fcn(x, p);
			if (j == lastRealIndexP)
				GG_Top[i][j] = Dx * wVal * boundaryVal;
			else if (wVal >= 0)
				GG_Top[i][j] = Dx * wVal * sl[i][j];
			else
				GG_Top[i][j] = Dx * wVal * sl[i][j + 1];

			/*
			 * FF fluxes
			 */
			p = getCellCenterP(i, j);

			// FF left
			x = getCellLeftX(i, j);
			uVal = u_fcn(x, p);
			if (i == 1) {
				FF_Left[i][j] = -Dp * u_fcn(x, p) * boundaryVal;
			}
			else if (uVal >= 0) {
				FF_Left[i][j] = -Dp * uVal * sl[i - 1][j];
			}
			else {
				FF_Left[i][j] = -Dp * uVal * sl[i][j];
			}

			// FF right
			x = getCellRightX(i, j);
			uVal = u_fcn(x, p);
			if (i == lastRealIndexX)
				FF_Right[i][j] = Dp * u_fcn(x, p) * boundaryVal;
			else if (uVal >= 0)
				FF_Right[i][j] = Dp * uVal * sl[i][j];
			else
				FF_Right[i][j] = Dp * uVal * sl[i + 1][j];
		}
}

void calcFluxes() {
	//return calcFluxes_Godunov_Dirichlet();
	return calcFluxes_ClassFV_Dirichlet();
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Time Steps
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void forwardEuler() {
	enforceInitCond();
	for (int tt = 0; tt < numTimeSteps; tt++) {
		double t = tt * Dt;
		//printf("\rProgress: %1.1f%%", t / finalTime * 100);
		calcFluxes();
		for (int i = 1; i < lastGhostIndexX; i++)
			for (int j = 1; j < lastGhostIndexP; j++) {
				double T = sl[i][j], x = getCellCenterX(i, j), p = getCellCenterP(i, j);
				double RHS = source_Test1(T, x, p, t) -
						(GG[i][j] - GG[i][j - 1] + FF[i][j] - FF[i - 1][j]) / cellVol;
				sl[i][j] += Dt * RHS;
			}
	}
	printf("\nCompleted: 100%%\n");
}

void forwardEuler2() {
	enforceInitCond();
	for (int tt = 0; tt < numTimeSteps; tt++) {
		double t = tt * Dt;
		//printf("\rProgress: %1.1f%%", t / finalTime * 100);
		calcFluxes();
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
