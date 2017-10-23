/*
 * QuadCell.h
 *
 *  Created on: Oct 23, 2017
 *      Author: chuckjia
 */

#ifndef QUADCELL_H_
#define QUADCELL_H_
#include "Projection.h"

double a1_quadCell_cache;
double a2_quadCell_cache[Nx + 1][Np + 1], a3_quadCell_cache[Nx + 1][Np + 1],
a4_quadCell_cache[Nx + 1][Np + 1];

void fillCache_quadCell() {
	a1_quadCell_cache = 0.25;
	double a1 = a1_quadCell_cache, one_minus_a1 = 1 - a1;
	for (int i = 0; i <= Nx; ++i) {
		double xThis = getCellCenterX(i), xRight = getCellRightX(i),
				xNext = getCellCenterX(i + 1);
		double d22 = xThis - xNext;
		double e2 = xRight - a1 * xThis - one_minus_a1 * xNext;
		for (int j = 1; j <= Np; ++j) {
			double pThis = getCellCenterP(i, j), pThisTopRight = getCellTopRightP(i, j),
					pUpperNext = getCellCenterP(i, j + 1),
					pRightNext = getCellCenterP(i + 1, j),
					pUpperRightNext = getCellCenterP(i + 1, j + 1);
			double d32 = pUpperNext - pRightNext, d33 = pUpperRightNext - pRightNext,
					e3 = pThisTopRight - a1 * pThis - one_minus_a1 * pRightNext;
			double a3 = e2 / d22,
					a4 = (e3 - d32 * a3) / d33,
					a2 = one_minus_a1 - a3 - a4;
			a2_quadCell_cache[i][j] = a2;
			a3_quadCell_cache[i][j] = a3;
			a4_quadCell_cache[i][j] = a4;
		}
	}
}

double e11_diagMatInv_quadCell_cache, e21_diagMatInv_quadCell_cache;
double e12_diagMatInv_quadCell_cache[Nx + 1][Np + 1],
e22_diagMatInv_quadCell_cache[Nx + 1][Np + 1];

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

double getCellTopRightU(int i, int j) {
	return a1_quadCell_cache * u[i][j] + a2_quadCell_cache[i][j] * u[i + 1][j] +
			a3_quadCell_cache[i][j] * u[i][j + 1] + a4_quadCell_cache[i][j] * u[i + 1][j + 1];
}

double getGradhU_x(int i, int j) {
	return e11_diagMatInv_quadCell_cache * (getCellTopRightU(i, j) - getCellTopRightU(i - 1, j)) +
			e12_diagMatInv_quadCell_cache[i][j] * (u[i][j + 1] - u[i][j]);
}

double getGradhU_p(int i, int j) {
	return e22_diagMatInv_quadCell_cache[i][j] * (u[i][j + 1] - u[i][j]);
}

double getCellTopRightT(int i, int j) {
	return a1_quadCell_cache * u[i][j] + a2_quadCell_cache[i][j] * u[i + 1][j] +
			a3_quadCell_cache[i][j] * u[i][j + 1] + a4_quadCell_cache[i][j] * u[i + 1][j + 1];
}

double getGradhT_x(int i, int j) {
	return e11_diagMatInv_quadCell_cache * (getCellTopRightT(i, j) - getCellTopRightT(i - 1, j)) +
			e12_diagMatInv_quadCell_cache[i][j] * (T[i][j + 1] - T[i][j]);
}

double getGradhT_p(int i, int j) {
	return e22_diagMatInv_quadCell_cache[i][j] * (T[i][j + 1] - T[i][j]);
}

void setQuadCells() {
	fillCache_quadCell();
	fillCache_diagMatInv_quadCell();
}

#endif /* QUADCELL_H_ */
