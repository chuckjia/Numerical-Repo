/*
 * QuadCell.h
 *
 *  Created on: Oct 23, 2017
 *      Author: chuckjia
 *
 *  This file contains global arrays and functions that computes quadrilateral cell interpolations.
 */

#ifndef E_QUADCELL_H_
#define E_QUADCELL_H_
#include "D_Mesh.h"

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Interpolation Coefficients for The Quadrilateral Cells
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Cache for the interpolation coefficients on The quadrilateral cells. See (3.13)
double a1_quad;
double a2_quad_[Nx + 1][Np + 1], a3_quad_[Nx + 1][Np + 1], a4_quad_[Nx + 1][Np + 1];

// Calculate the coefficients a1 - a4. We set a1 to be 0.25.
void fillCache_quadCell() {
	a1_quad = 0.25;
	double a1 = a1_quad;
	for (int i = 0; i <= Nx; ++i) {
		double thisCell_centerX = getCellCenterX(i), thisCell_rightX = getCellRightX(i), rightCell_centerX = getCellCenterX(i + 1);
		double d22 = thisCell_centerX - rightCell_centerX,
				e2 = thisCell_rightX - a1 * thisCell_centerX - (1 - a1) * rightCell_centerX;
		for (int j = 1; j <= Np; ++j) {
			double thisCell_centerP = getCellCenterP(i, j), thisCell_topRightP = getCellTopRightP(i, j),
					upperCell_centerP = getCellCenterP(i, j + 1),
					rightCell_centerP = getCellCenterP(i + 1, j),
					upperRightCell_centerP = getCellCenterP(i + 1, j + 1);
			double d32 = upperCell_centerP - rightCell_centerP,
					d33 = upperRightCell_centerP - rightCell_centerP,
					e3 = thisCell_topRightP - a1 * thisCell_centerP - (1 - a1) * rightCell_centerP;
			double a3 = e2 / d22,
					a4 = (e3 - d32 * a3) / d33,
					a2 = 1 - a1 - a3 - a4;
			a2_quad_[i][j] = a2;
			a3_quad_[i][j] = a3;
			a4_quad_[i][j] = a4;
		}
	}
}

// Cache for the matrix M_{i,j+1/2}, which represents the diagonal of cell(i, j+1/2).
// Only the inverse of the matrix is actually stored
double e11_MInv_quad, e21_MInv_quad;
double e12_MInv_quad_[Nx + 1][Np + 1], e22_MInv_quad_[Nx + 1][Np + 1];

// Calculate the matrices for quadrilateral cell diagonals. Only inverses are stored in cache
void fillCache_diagMatInv_quadCell() {
	e11_MInv_quad = 1. / Dx;
	e21_MInv_quad = 0;
	for (int i = 0; i <= Nx; ++i)
		for (int j = 0; j <= Np; ++j) {  // Variables a, b, and d correspond to the variables used in the notes
			double b = getCellTopRightP(i, j) - getCellTopLeftP(i, j),
					d = getCellCenterP(i, j + 1) - getCellCenterP(i, j);
			e12_MInv_quad_[i][j] = -b / (Dx * d);
			e22_MInv_quad_[i][j] = 1. / d;
		}
}

// Return u_{i+1/2, j+1/2}, i.e. the interpolated value of u in the quadrilateral cells. See (3.14)
double getCellTopRightU(int i, int j) {
	return a1_quad * u_[i][j] + a2_quad_[i][j] * u_[i + 1][j] +
			a3_quad_[i][j] * u_[i][j + 1] + a4_quad_[i][j] * u_[i + 1][j + 1];
}

/*
 * Getters for gradient_h u_h over the quadrilateral cell (i,j+1/2). See (3.15)
 *   - The two components for the results are separated because only getGradhU_x will be used
 */

double getGradhU_x(int i, int j) {
	return e11_MInv_quad * (getCellTopRightU(i, j) - getCellTopRightU(i - 1, j)) +
			e12_MInv_quad_[i][j] * (u_[i][j + 1] - u_[i][j]);
}

double getGradhU_p(int i, int j) {
	return e22_MInv_quad_[i][j] * (u_[i][j + 1] - u_[i][j]);  // Assuming the (2, 1) element of the inverse matrix is 0
}

// Return T_{i+1/2, j+1/2}, i.e. the interpolated value of T in the quadrilateral cells
double getCellTopRightT(int i, int j) {
	return a1_quad * T_[i][j] + a2_quad_[i][j] * T_[i + 1][j] +
			a3_quad_[i][j] * T_[i][j + 1] + a4_quad_[i][j] * T_[i + 1][j + 1];
}

/*
 * Getters for gradient_h T_h over the quadrilateral cell (i,j+1/2)
 *   - The two components for the results are separated because only getGradhT_x will be used
 */

double getGradhT_x(int i, int j) {
	return e11_MInv_quad * (getCellTopRightT(i, j) - getCellTopRightT(i - 1, j)) +
			e12_MInv_quad_[i][j] * (T_[i][j + 1] - T_[i][j]);
}

double getGradhT_p(int i, int j) {
	return e22_MInv_quad_[i][j] * (T_[i][j + 1] - T_[i][j]);
}


/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Set All Parameters And Calculate All Cache Values For the Quadrilateral Cells
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void setQuadCell() {
	fillCache_quadCell();
	fillCache_diagMatInv_quadCell();
}

#endif /* E_QUADCELL_H_ */
