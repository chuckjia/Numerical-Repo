/*
 * Fluxes.h
 *
 *  Created on: Jun 10, 2017
 *      Author: chuckjia
 */

#ifndef FLUXES_H_
#define FLUXES_H_
#include "ModelConditions.h"

// Initialize fluxes
double Hx[numCellsXDir][numCellsPDir][2];
double Hp[numCellsXDir][numCellsPDir][2];

// Initialize the reconstructions
double uX[numCellsXDir][numCellsPDir][2];
double uP[numCellsXDir][numCellsPDir][2];

void reconstrSoln() {
	for (int j = 1; j < numGridPtsXDir; j++) {
		double DpVal = getDp(j);
		for (int k = 1; k < numGridPtsPDir; k++) {
			double solnThis1 = soln[j][k][0], solnThis2 = soln[j][k][1],
					solnRight1 = soln[j + 1][k][0], solnRight2 = soln[j + 1][k][1],
					solnLeft1 = soln[j - 1][k][0], solnLeft2 = soln[j - 1][k][1],
					solnUpper1 = soln[j][k + 1][0], solnUpper2 = soln[j][k + 1][1],
					solnLower1 = soln[j][k - 1][0], solnLower2 = soln[j][k - 1][1];
			uX[j][k][0] = minmod3(theta_CONST * (solnRight1 - solnThis1) / Dx,
					0.5 * (solnRight1 - solnLeft1) / Dx,
					theta_CONST * (solnThis1 - solnLeft1) / Dx);
			uX[j][k][1] = minmod3(theta_CONST * (solnRight2 - solnThis2) / Dx,
					0.5 * (solnRight2 - solnLeft2) / Dx,
					theta_CONST * (solnThis2 - solnLeft2) / Dx);
			uP[j][k][0] = minmod3(theta_CONST * (solnUpper1 - solnThis1) / DpVal,
					0.5 * (solnUpper1 - solnLower1) / DpVal,
					theta_CONST * (solnThis1 - solnLower1) / DpVal);
			uP[j][k][1] = minmod3(theta_CONST * (solnUpper2 - solnThis2) / DpVal,
					0.5 * (solnUpper2 - solnLower2) / DpVal,
					theta_CONST * (solnThis2 - solnLower2) / DpVal);
		}
	}
	// When j = 0 or numCellsPDir - 1
	int indRange[2] = {0, numCellsPDir - 1};
	for (int ii = 0; ii < 2; ii++) {
		int jVal = indRange[ii];
		for (int k = 1; k < numGridPtsPDir; k++) {
			uX[jVal][k][0] = 0;
			uX[jVal][k][1] = 0;
			double DpVal = getDp(jVal);
			double solnThis1 = soln[jVal][k][0], solnThis2 = soln[jVal][k][1],
					solnUpper1 = soln[jVal][k + 1][0], solnUpper2 = soln[jVal][k + 1][1],
					solnLower1 = soln[jVal][k - 1][0], solnLower2 = soln[jVal][k - 1][1];
			uP[jVal][k][0] = minmod3(theta_CONST * (solnUpper1 - solnThis1) / DpVal,
					0.5 * (solnUpper1 - solnLower1) / DpVal,
					theta_CONST * (solnThis1 - solnLower1) / DpVal);
			uP[jVal][k][1] = minmod3(theta_CONST * (solnUpper2 - solnThis2) / DpVal,
					0.5 * (solnUpper2 - solnLower2) / DpVal,
					theta_CONST * (solnThis2 - solnLower2) / DpVal);
		}
	}
	// When i = 0 or numCellsXDir - 1
	indRange[1] = numCellsXDir - 1;
	for (int ii = 0; ii < 2; ii++) {
		int kVal = indRange[ii];
		for (int j = 1; j < numGridPtsPDir; j++) {
			double solnThis1 = soln[j][kVal][0], solnThis2 = soln[j][kVal][1],
					solnRight1 = soln[j + 1][kVal][0], solnRight2 = soln[j + 1][kVal][1],
					solnLeft1 = soln[j - 1][kVal][0], solnLeft2 = soln[j - 1][kVal][1];
			uX[j][kVal][0] = minmod3(theta_CONST * (solnRight1 - solnThis1) / Dx,
					0.5 * (solnRight1 - solnLeft1) / Dx,
					theta_CONST * (solnThis1 - solnLeft1) / Dx);
			uX[j][kVal][1] = minmod3(theta_CONST * (solnRight2 - solnThis2) / Dx,
					0.5 * (solnRight2 - solnLeft2) / Dx,
					theta_CONST * (solnThis2 - solnLeft2) / Dx);
			uP[j][kVal][0] = 0;
			uP[j][kVal][1] = 0;
		}
	}
	// Corners
	uX[0][0][0] = 0; uX[0][0][1] = 0; uP[0][0][0] = 0; uP[0][0][1] = 0;
	int jTemp = numCellsXDir - 1, kTemp = numCellsPDir - 1;
	uX[jTemp][0][0] = 0; uX[jTemp][0][1] = 0; uP[jTemp][0][0] = 0; uP[jTemp][0][1] = 0;
	uX[0][kTemp][0] = 0; uX[0][kTemp][1] = 0; uP[0][kTemp][0] = 0; uP[0][kTemp][1] = 0;
	uX[jTemp][kTemp][0] = 0; uX[jTemp][kTemp][1] = 0;
	uP[jTemp][kTemp][0] = 0; uP[jTemp][kTemp][1] = 0;
}

void reconstrFcn(double res[2], int j, int k, double x, double p) {
	double xCenter = getCenterX(j, k), pCenter = getCenterP(j, k);
	double uThis1 = soln[j][k][0], uThis2 = soln[j][k][1];
	double xTerm = x - xCenter; // The term x - x_j
	double pTerm = p - pCenter; // The term y - y_k
	res[0] = uThis1 + uX[j][k][0] * xTerm + uP[j][k][0] * pTerm;
	res[1] = uThis2 + uX[j][k][1] * xTerm + uP[j][k][1] * pTerm;
}

void calcFluxes() {
	reconstrSoln();
	for (int j = 0; j < numCellsXDir; j++) {
		double DpVal = getDp(j);
		for (int k = 0; k < numCellsPDir; k++) {
			// The center of this cell
			double xCenter = getCenterX(j, k), pCenter = getCenterP(j, k);

			// Reconstructed values at the interfaces
			double solnN[2], solnSNext[2], solnE[2], solnWNext[2];
			double pTemp = 0.5 * (getTopLeftP(j, k) + getTopRightP(j, k));
			reconstrFcn(solnN, j, k, xCenter, pTemp);
			reconstrFcn(solnSNext, j, k + 1, xCenter, pTemp);
			double xTemp = getCellRightX(j, k);
			reconstrFcn(solnE, j, k, xTemp, pCenter);
			reconstrFcn(solnWNext, j + 1, k, xTemp, pCenter);

			// The speeds of the interface
			// Speeds a+ and a-
			double eval1 = uFcn(solnWNext[0], solnWNext[1]),
					eval2 = uFcn(solnE[0], solnE[1]); // Eigenvalues
			double aPlus = max3(eval1, eval2, 0);
			double aMinus = min3(eval1, eval2, 0);
			// Speeds b+ and b-, 1st component
			eval1 = omegaFcn(solnSNext[0], solnSNext[1]);
			eval2 = omegaFcn(solnN[0], solnN[1]);
			double bPlus = max3(eval1, eval2, 0);
			double bMinus = min3(eval1, eval2, 0);

			// Calculate the numerical fluxes: Hx
			double plusVal = aPlus, minusVal = aMinus;
			xTemp = getCellRightX(j, k), pTemp = pCenter;
			// Here solnVal = u_{j, k}^E, solnNextVal = u_{j + 1,k}^W
			// 1st component of Hx
			double solnVal = solnE[0], solnNextVal = solnWNext[0];
			Hx[j][k][0] = (plusVal * fFcn(solnVal, xTemp, pTemp) -
					minusVal * fFcn(solnNextVal, xTemp, pTemp) +
					plusVal * minusVal * (solnNextVal - solnVal))
					/ (plusVal - minusVal);
			// 2nd component of Hx
			solnVal = solnE[1], solnNextVal = solnWNext[1];
			Hx[j][k][1] = (plusVal * fFcn(solnVal, xTemp, pTemp) -
					minusVal * fFcn(solnNextVal, xTemp, pTemp) +
					plusVal * minusVal * (solnNextVal - solnVal))
					/ (plusVal - minusVal);
			// Calculate the numerical fluxes: Hp
			plusVal = bPlus, minusVal = bMinus;
			xTemp = xCenter, pTemp = 0.5 * (getTopLeftP(j, k) + getTopRightP(j, k));
			// Here solnVal = u_{j, k}^N, solnNextVal = u_{j + 1,k}^S
			// 1st component of Hp
			solnVal = solnN[0], solnNextVal = solnSNext[0];
			Hp[j][k][0] = (plusVal * fFcn(solnVal, xTemp, pTemp) -
					minusVal * fFcn(solnNextVal, xTemp, pTemp) +
					plusVal * minusVal * (solnNextVal - solnVal))
					/ (plusVal - minusVal);
			// 2nd component of Hp
			solnVal = solnN[1], solnNextVal = solnSNext[1];
			Hp[j][k][1] = (plusVal * fFcn(solnVal, xTemp, pTemp) -
					minusVal * fFcn(solnNextVal, xTemp, pTemp) +
					plusVal * minusVal * (solnNextVal - solnVal))
					/ (plusVal - minusVal);
		}
	}
}

#endif /* FLUXES_H_ */
