/*
 * w.h
 *
 *  Created on: Oct 21, 2017
 *      Author: chuckjia
 */

#ifndef W_H_
#define W_H_
#include "Projection.h"

double a1_quadCell;
double a2_quadCell[Nx + 1][Np + 1], a3_quadCell[Nx + 1][Np + 1], a4_quadCell[Nx + 1][Np + 1];

void fillCache_quadCell() {
	a1_quadCell = 0.25;
	double a1 = a1_quadCell, one_minus_a1 = 1 - a1;
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
			a2_quadCell[i][j] = a2;
			a3_quadCell[i][j] = a3;
			a4_quadCell[i][j] = a4;
		}
	}
}

double getCellTopRight_u(int i, int j) {
	return a1_quadCell * u[i][j] + a2_quadCell[i][j] * u[i + 1][j] +
			a3_quadCell[i][j] * u[i][j + 1] + a4_quadCell[i][j] * u[i + 1][j + 1];
}

double e11_diagMatInv_quadCell = DxInv, e21_diagMatInv_quadCell = 0;
double e12_diagMatInv_quadCell[Nx + 1][Np + 1], e22_diagMatInv_quadCell[Nx + 1][Np + 1];

void fillCache_diagMatInv_quadCell() {
	for (int i = 0; i <= Nx; ++i)
		for (int j = 1; j <= Np; ++j) {
			double b = getCellTopRightP(i, j) - getCellTopLeftP(i, j),
					dInv = 1 / (getCellCenterP(i, j + 1) - getCellCenterP(i, j));
			e12_diagMatInv_quadCell[i][j] = - b * DxInv * dInv;
			e22_diagMatInv_quadCell[i][j] = dInv;
		}
}

double get_uGradX(int i, int j) {
	return e11_diagMatInv_quadCell * (getCellTopRight_u(i, j) - getCellTopRight_u(i - 1, j)) +
			e12_diagMatInv_quadCell[i][j] * (u[i][j + 1] - u[i][j]);
}

double get_uGradP(int i, int j) {
	return e22_diagMatInv_quadCell[i][j] * (u[i][j + 1] - u[i][j]);
}

void calc_w() {
	for (int i = 1; i <= Nx; ++i)
		for (int j = 0; j < Np; ++j)
			w[i][j + 1] = w[i][j] - (getCellCenterP(i, j + 1) - getCellCenterP(i, j)) *
					get_uGradX(i, j);

}

void set_w() {
	fillCache_quadCell();
	fillCache_diagMatInv_quadCell();
}

#endif /* W_H_ */
