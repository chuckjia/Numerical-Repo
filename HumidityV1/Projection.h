/*
 * Projection.h
 *
 *  Created on: Oct 19, 2017
 *      Author: chuckjia
 */

#ifndef PROJECTION_H_
#define PROJECTION_H_
#include "QuadCell.h"

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Projection Method For u
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Coefficients for the linear system (3.33) + (3.35) in the projection method. To get the
// corresponding matrix, we choose to put (3.35) as the last row. Then we directly reduce
// the last row to make the matrix upper triangular, and we solve the system backwards.

// Our a and b are different from those in the article:
// a_this[i][j] = a_orig[i][j] - b_orig[i][j] / Dx, b_this[i][j] = b_orig[i][j] / Dx.
// In our program here, a and b are chosen to be exactly the matrix entries in (3.33).
// Note: we only store the inverse of a.

// d values represent the last row of the matrix, corresponding to (3.35). They store the
// original last matrix row and the intermediate values in the calculation of a and b. More
// importantly, d stores the linear transformation on c to directly make the matrix upper
// triangular.

// During the life of this program, the cache values of a, b, and d will remain the same once
// calculated.

// Refer to notes for details.

double aInv_proj[Nx], b_proj_[Nx];  // Only values from 1 to Nx-1 are used
double d_proj_[Nx + 1]; // Only values from 1 to Nx are used
double lambdax_proj_[Nx + 1];  // Only values from 1 to Nx are used

// Calculate a and b values
void fillCache_ab_proj() {
	for (int i = 1; i < Nx; ++i) {
		double x = getCellCenterX(i);
		double bThis = ((*pB_fptr)(x) - pA) * DxInv;
		b_proj_[i] = bThis;
		aInv_proj[i] = 1 / ((*pBxDer_fptr)(x) - bThis);
	}
}

// Calculate d values. It is separated from the calculation of a and b for testing reasons
void fillCache_d_proj() {
	d_proj_[1] = 1;
	for (int i = 1; i < Nx; ++i) {
		double temp = d_proj_[i] * aInv_proj[i];
		d_proj_[i] = temp;
		d_proj_[i + 1] = 1 - b_proj_[i] * temp;
	}
}

// Calculate the c_i's as in (3.34).
// In this program, we do not have an array for c_i. Instead, we use lambda_x[] to hold the c_i
// values. Later after applying Gaussian elimination, lambda_x[] will naturally hold the
// lambda_x values as the solution of the linear system.
void calc_c_proj() {
	// Temporarily store the integer int_pA^pB(x_i)\tilde{u}dp in c_i. The last c_Nx will not
	// be updated to 0, but since we directly assign value to c_Nx in the Gaussian elimination,
	// this will not affect the calculation.
	for (int i = 1; i <= Nx; ++i) {
		double sum = 0;
		for (int j = 1; j <= lastRealIndexP; j++)
			sum += u_[i][j];
		sum *= getCellCenterDp(i);
		lambdax_proj_[i] = sum;
	}
	// Now calculate the real c_proj values
	for (int i = 1; i < Nx; ++i)
		lambdax_proj_[i] = (lambdax_proj_[i + 1] - lambdax_proj_[i]) * DxInv;
}

// Perform Gaussian elimination to calculate lambda_x
void calcLambdax_proj() {
	calc_c_proj();
	double sum = 0;
	for (int i = 1; i < Nx; ++i)
		sum += lambdax_proj_[i] * d_proj_[i];
	lambdax_proj_[Nx] = -sum / d_proj_[Nx];
	for (int i = Nx - 1; i >= 1; --i)
		lambdax_proj_[i] = (lambdax_proj_[i] - b_proj_[i] * lambdax_proj_[i + 1]) * aInv_proj[i];
}

// Perform the projection method on uTilde to calculate u
void (*projU_fptr)();

// Perform the projection method on uTilde to calculate u
void projU_orig() {
	calcLambdax_proj();
	for (int i = 1; i < Nx; ++i) {
		double lambda_x = lambdax_proj_[i];
		for (int j = 1; j < lastRealIndexP; ++j)
			u_[i][j] -= lambda_x;
	}
}

// For testing
void print_lambdax() {
	printf("\n");
	for (int i = 1; i <= Nx; ++i)
		printf("lambda_x[%d] = %1.2e  ", i, lambdax_proj_[i]);
	printf("\n");
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Set Projection Parameters and Fill Cache Values
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void setProjection() {
	fillCache_ab_proj();
	fillCache_d_proj();
	projU_fptr = &projU_orig;
}


#endif /* PROJECTION_H_ */
