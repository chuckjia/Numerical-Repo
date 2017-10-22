/*
 * Mesh.h
 *
 *  Created on: Sep 18, 2017
 *      Author: chuckjia
 */

#ifndef MESH_H_
#define MESH_H_
#include "Constants.h"

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Getter Functions For the Mesh
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Return the x coordinate of the left side of cell (i, j)
double getCellLeftX(int i, int j) {
	return x0 + (i - 1) * Dx;
}

// Return the x coordinate of the left side of cell (i, j) using only i index
double getCellLeftX(int i) {
	return x0 + (i - 1) * Dx;
}

// Return the x coordinate of the right side of cell (i, j)
double getCellRightX(int i, int j) {
	return x0 + i * Dx;
}

// Return the x coordinate of the right side of cell (i, j) using only i index
double getCellRightX(int i) {
	return x0 + i * Dx;
}

// Store the values of Dp. Its ith element is the Dp value at x_{i+1/2}
double DpValMat[numGridPtsX];

// Calculate all Dp values and store in DpValMat
double calcDpVal() {
	for (int i = 1; i < numGridPtsX; i++) {
		double x = getCellLeftX(i);
		DpValMat[i] = (pB(x) - pA) / Np;
	}
	// On ghost cell boundaries that are outside the domain
	DpValMat[0] = DpValMat[1];
	DpValMat[numGridPtsX - 1] = DpValMat[numGridPtsX - 2];
}

// Return the Dp value on the left side of cell (i, j)
double getCellLeftDpVal(int i, int j) {
	return DpValMat[i];
}

// Return the Dp value on the left side of cell (i, j), using only index i
double getCellLeftDpVal(int i) {
	return DpValMat[i];
}

// Return the Dp value on the right side of cell (i, j)
double getCellRightDpVal(int i, int j) {
	return DpValMat[i + 1];
}

// Return the Dp value on the left side of cell (i, j), using only index i
double getCellRightDpVal(int i) {
	return DpValMat[i + 1];
}

// Return the p coordinate of the bottom side of cell (i, j)
double getCellBottLeftP(int i, int j) {
	return pA + (j - 1) * getCellLeftDpVal(i);
}

// Return the p coordinate of the top side of cell (i, j)
double getCellTopLeftP(int i, int j) {
	return pA + j * getCellLeftDpVal(i);
}

// Return the p coordinate of the top side of cell (i, j)
double getCellBottRightP(int i, int j) {
	return pA + (j - 1) * getCellRightDpVal(i);
}

// Return the p coordinate of the top side of cell (i, j)
double getCellTopRightP(int i, int j) {
	return pA + j * getCellRightDpVal(i);
}

// Store the coordinates of cell centers
double cellCentersX[numCellsX];
double cellCentersP[numCellsX][numCellsP];

// Calculate cell centers and store them in cellCentersX and cellCentersP
void calcCellCenters() {
	for (int i = 0; i < numCellsX; i++) {
		double a = getCellLeftDpVal(i), b = getCellRightDpVal(i);
		cellCentersX[i] = x0 + i * Dx - Dx / 3 * (2 * a + b) / (a + b);
		for (int j = 0; j < numCellsP; j++) {
			double bottRightP = getCellBottRightP(i, j),
					c = bottRightP - getCellBottLeftP(i, j);
			cellCentersP[i][j] = bottRightP +
					(2 * a * c + a * a + c * b + a * b + b * b) / (a + b) / 3;
		}
	}
}

// Return the x coordinate of the center of cell (i, j)
double getCellCenterX(int i, int j) {
	return cellCentersX[i];
}

// Return the x coordinate of the center of cell (i, j), using only index i
double getCellCenterX(int i) {
	return cellCentersX[i];
}

// Return the p coordinate of the center of cell (i, j)
double getCellCenterP(int i, int j) {
	return cellCentersP[i][j];
}

// Store the volume of all cells
double cellVol[numCellsX];


// Calculate all cell volume and store them in cellVol
void calcCellVol() {
	for (int i = 0; i < numCellsX; i++) {
		double a = getCellLeftDpVal(i), b = getCellRightDpVal(i);
		cellVol[i] = (a + b) * Dx / 2;
	}
}

// Return the volume of cell (i, j)
double getCellVol(int i, int j) {
	return cellVol[i];
}

// Return the volume of cell (i, j), using only index i
double getCellVol(int i) {
	return cellVol[i];
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Model Functions
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Coefficients c1 and c2 used in the u function
const double c1_uFcn = M_PI / pB, c2_uFcn = 2 * M_PI / xf;

// Velocity function u
double u_fcn(double x, double p) {
	return 7.5 + cos(c1_uFcn * p) * cos(c2_uFcn * x);
}

// x derivative of the u function
double u_xDer_fcn(double x, double p) {
	return -cos(c1_uFcn * p) * c2_uFcn * sin(c2_uFcn * x);
}

// Coefficients c1 and c2 used in the w function
const double c1_wFcn = M_PI / pB, c2_wFcn = 2 * M_PI / xf;

// Velocity function w
double w_fcn(double x, double p) {
	return sin(c1_wFcn * p) * cos(c2_wFcn * x);
}

// p derivative of the w function
double w_pDer_fcn(double x, double p) {
	return c1_wFcn * cos(c1_wFcn * p) * cos(c2_wFcn * x);
}

#endif /* MESH_H_ */
