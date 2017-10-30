/*
 * QuadCell.h
 *
 *  Created on: Oct 23, 2017
 *      Author: chuckjia
 */

#ifndef QUADCELL_H_
#define QUADCELL_H_
#include "Projection.h"

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Interpolation Coefficients For The Quadrilateral Cells
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Cache for the interpolation coefficients on The quadrilateral cells. See (3.13).
double a1_quadCell_cache;
double a2_quadCell_cache[Nx + 1][Np + 1], a3_quadCell_cache[Nx + 1][Np + 1],
a4_quadCell_cache[Nx + 1][Np + 1];

// Calculate the coefficients a1 - a4. We set a1 to be 0.25.
void fillCache_quadCell() {
	a1_quadCell_cache = 0.25;
	double a1 = a1_quadCell_cache, one_minus_a1 = 1 - a1;
	for (int i = 0; i <= Nx; ++i) {
		double xThisCenter = getCellCenterX(i), xCellRight = getCellRightX(i),
				xNextCenter = getCellCenterX(i + 1);
		double d22 = xThisCenter - xNextCenter,
				e2 = xCellRight - a1 * xThisCenter - one_minus_a1 * xNextCenter;
		for (int j = 1; j <= Np; ++j) {
			double pThisCenter = getCellCenterP(i, j), pThisTopRight = getCellTopRightP(i, j),
					pUpperNextCenter = getCellCenterP(i, j + 1),
					pRightNextCenter = getCellCenterP(i + 1, j),
					pUpperRightNextCenter = getCellCenterP(i + 1, j + 1);
			double d32 = pUpperNextCenter - pRightNextCenter,
					d33 = pUpperRightNextCenter - pRightNextCenter,
					e3 = pThisTopRight - a1 * pThisCenter - one_minus_a1 * pRightNextCenter;
			double a3 = e2 / d22,
					a4 = (e3 - d32 * a3) / d33,
					a2 = one_minus_a1 - a3 - a4;
			a2_quadCell_cache[i][j] = a2;
			a3_quadCell_cache[i][j] = a3;
			a4_quadCell_cache[i][j] = a4;
		}
	}
}

// Cache for the matrix M_{i,j+1/2}, which represents the diagonal of cell(i, j+1/2).
// Only the inverse of the matrix is actually stored
double e11_diagMatInv_quadCell_cache, e21_diagMatInv_quadCell_cache;
double e12_diagMatInv_quadCell_cache[Nx + 1][Np + 1],
e22_diagMatInv_quadCell_cache[Nx + 1][Np + 1];

// Calculate the matrices for quadrilateral cell diagonals. Only inverses are stored in cache
void fillCache_diagMatInv_quadCell() {
	e11_diagMatInv_quadCell_cache = DxInv;
	e21_diagMatInv_quadCell_cache = 0;
	for (int i = 0; i <= Nx; ++i)
		for (int j = 1; j <= Np; ++j) {
			double b = getCellTopRightP(i, j) - getCellTopLeftP(i, j),
					dInv = 1 / (getCellCenterP(i, j + 1) - getCellCenterP(i, j));
			e12_diagMatInv_quadCell_cache[i][j] = - b * DxInv * dInv;
			e22_diagMatInv_quadCell_cache[i][j] = dInv;
		}
}

// Return the interpolated value of u in the quadrilateral cells
double getCellTopRightU(int i, int j) {
	return a1_quadCell_cache * u_sl[i][j] + a2_quadCell_cache[i][j] * u_sl[i + 1][j] +
			a3_quadCell_cache[i][j] * u_sl[i][j + 1] + a4_quadCell_cache[i][j] * u_sl[i + 1][j + 1];
}

// Getters for gradient_h u_h over the quadrilateral cell (i,j+1/2)
// The two components for the results are separated because only getGradhU_x will be used

double getGradhU_x(int i, int j) {
	return e11_diagMatInv_quadCell_cache * (getCellTopRightU(i, j) - getCellTopRightU(i - 1, j)) +
			e12_diagMatInv_quadCell_cache[i][j] * (u_sl[i][j + 1] - u_sl[i][j]);
}

// This getter assumes that the (2, 1) element of the inverse matrix is 0
double getGradhU_p(int i, int j) {
	return e22_diagMatInv_quadCell_cache[i][j] * (u_sl[i][j + 1] - u_sl[i][j]);
}

// Return the interpolated value of T in the quadrilateral cells
double getCellTopRightT(int i, int j) {
	return a1_quadCell_cache * T_sl[i][j] + a2_quadCell_cache[i][j] * T_sl[i + 1][j] +
			a3_quadCell_cache[i][j] * T_sl[i][j + 1] + a4_quadCell_cache[i][j] * T_sl[i + 1][j + 1];
}

// Getters for gradient_h T_h over the quadrilateral cell (i,j+1/2)
// The two components for the results are separated because only getGradhT_x will be used

double getGradhT_x(int i, int j) {
	return e11_diagMatInv_quadCell_cache * (getCellTopRightT(i, j) - getCellTopRightT(i - 1, j)) +
			e12_diagMatInv_quadCell_cache[i][j] * (T_sl[i][j + 1] - T_sl[i][j]);
}

// Getters for gradient_h Ts_h over the quadrilateral cell (i,j+1/2)

double getGradhT_p(int i, int j) {
	return e22_diagMatInv_quadCell_cache[i][j] * (T_sl[i][j + 1] - T_sl[i][j]);
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Set All Parameters And Calculate All Cache Values For the Quadrilateral Cells
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void setQuadCells() {
	fillCache_quadCell();
	fillCache_diagMatInv_quadCell();
}

#endif /* QUADCELL_H_ */
