/*
 * Projection.h
 *
 *  Created on: Oct 19, 2017
 *      Author: chuckjia
 */

#ifndef G_PROJECTION_H_
#define G_PROJECTION_H_
#include "F_QuadCell.h"

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

double a_proj[Nx + 2], b_proj[Nx + 2];  // Only values from 1 to Nx-1 are used
double d_proj[Nx + 2]; // Only values from 1 to Nx are used
double lambda_x_proj[Nx + 2];  // Only values from 1 to Nx are used

int first_proj = 0, last_proj = Nx;

// Calculate a and b values
void fillCache_ab_proj() {
	for (int i = first_proj; i < last_proj; ++i) {  // Last index is last_proj - 1, see (3.33)
		double x = getCellCenterX(i);
		double bThis = (*pB_fptr)(x) - pA;
		b_proj[i] = bThis;
		a_proj[i] = (*pBxDer_fptr)(x) * Dx - bThis;
	}
}

// Calculate d values. It is separated from the calculation of a and b for testing reasons
void fillCache_d_proj() {
	d_proj[first_proj] = 1;
	for (int i = first_proj; i < last_proj; ++i) {  // Only up to (last index - 1)
		double temp = d_proj[i] / a_proj[i];
		d_proj[i] = temp;
		d_proj[i + 1] = 1 - b_proj[i] * temp;
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
	for (int i = first_proj; i <= last_proj; ++i) {  // Largest index is last_proj, see (3.33). In c_i, we need integral of u at x_{i+1}
		double sum = 0;
		for (int j = 1; j <= Np; j++)
			sum += u_[i][j];
		sum *= getCellCenterDp(i);
		lambda_x_proj[i] = sum;
	}
	// Now calculate the real c_proj values
	for (int i = first_proj; i < last_proj; ++i)  // Only up to (last index - 1)
		lambda_x_proj[i] = lambda_x_proj[i + 1] - lambda_x_proj[i];
}

// Perform Gaussian elimination to calculate lambda_x
void calcLambdax_gaussElim_proj() {
	calc_c_proj();
	double sum = 0;
	for (int i = first_proj; i < last_proj; ++i)
		sum += lambda_x_proj[i] * d_proj[i];
	lambda_x_proj[last_proj] = -sum / d_proj[last_proj];
	for (int i = last_proj - 1; i >= first_proj; --i)
		lambda_x_proj[i] = (lambda_x_proj[i] - b_proj[i] * lambda_x_proj[i + 1]) / a_proj[i];
}

// Perform the projection method on uTilde to calculate u
void (*projU_fcnPtr)();

// Perform the projection method on uTilde to calculate u
void projU_orig() {
	calcLambdax_gaussElim_proj();
	for (int i = first_proj; i <= last_proj; ++i) {
		double lambda_x = lambda_x_proj[i];
		for (int j = 1; j <= Np; ++j)
			u_[i][j] -= lambda_x;
	}
}

// For testing
void print_lambdax() {
	printf("\n");
	for (int i = first_proj; i <= Nx; ++i)
		printf("lambda_x[%d] = %1.2e  ", i, lambda_x_proj[i]);
	printf("\n");
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Set Projection Parameters and Fill Cache Values
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void setProjection() {
	fillCache_ab_proj();
	fillCache_d_proj();
	projU_fcnPtr = &projU_orig;
}


#endif /* G_PROJECTION_H_ */
