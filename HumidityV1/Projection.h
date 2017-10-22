/*
 * Projection.h
 *
 *  Created on: Oct 19, 2017
 *      Author: chuckjia
 */

#ifndef PROJECTION_H_
#define PROJECTION_H_
#include "Godunov.h"

double aInv_proj_cache[Nx], b_proj_cache[Nx];  // Only values from 1 to Nx-1 are used
double d_proj_cache[Nx + 1]; // Only values from 1 to Nx are used
double lambda_x[Nx + 1];  // Only values from 1 to Nx are used

void fillCache_ab_proj() {
	for (int i = 1; i < Nx; ++i) {
		double x = getCellCenterX(i);
		double temp = (*pB_fcnPtr)(x) * DxInv;
		b_proj_cache[i] = temp;
		aInv_proj_cache[i] = 1 / ((*pBxDer_fcnPtr)(x) - temp);
	}
}

// Calculate d values, separated for testing reasons
void fillCache_d_proj() {
	d_proj_cache[1] = 1;
	for (int i = 1; i < Nx; ++i) {
		double temp = d_proj_cache[i] * aInv_proj_cache[i];
		d_proj_cache[i] = temp;
		d_proj_cache[i + 1] = 1 - b_proj_cache[i] * temp;
	}
}

// Calculate the c_i's for the Gaussian elimination. For notations, see (3.33) and (3.34). In
// this program, lambda_x[] holds the c values at first. After applying Gaussian elimination, the
// lambda_x[] will hold the real lambda_x values.
void calc_c_proj() {
	// Temporarily store the integer int_pA^pB(x_i)\tilde{u}dp in c_i. The last c_Nx will not
	// be updated to 0, but since we directly assign value to c_Nx in the Gaussian elimination,
	// this is not a concern.
	for (int i = 1; i <= Nx; ++i) {
		double sum = 0;
		double x = getCellCenterX(i);
		for (int j = 1; j <= lastRealIndexP; j++)
			sum += u[i][j];
		sum *= getCellCenterDp(i);
		lambda_x[i] = sum;
	}
	// Now calculate the real c_proj values
	for (int i = 1; i < Nx; ++i)
		lambda_x[i] = (lambda_x[i + 1] - lambda_x[i]) * DxInv;
}

void gaussElim_proj() {
	calc_c_proj();
	double sum = 0;
	for (int i = 1; i < Nx; ++i)
		sum += lambda_x[i] * d_proj_cache[i];
	lambda_x[Nx] = -sum / d_proj_cache[Nx];
	for (int i = Nx - 1; i >= 1; --i)
		lambda_x[i] = (lambda_x[i] - b_proj_cache[i] * lambda_x[i + 1]) * aInv_proj_cache[i];
}

void setProjection() {
	fillCache_ab_proj();
	fillCache_d_proj();
}


#endif /* PROJECTION_H_ */
