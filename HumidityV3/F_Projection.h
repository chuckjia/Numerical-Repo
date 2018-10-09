/*
 * Projection.h
 *
 *  Created on: Oct 19, 2017
 *      Author: chuckjia
 *
 *  Global arrays and functions for the projection method.
 */

#ifndef F_PROJECTION_H_
#define F_PROJECTION_H_
#include "E_QuadCell.h"

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Projection Method For u
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void (*projU_fptr)();  // Function pointer for the projection method

int first_i_proj = 1, last_i_proj = Nx;  // Range of the index for the projection method

// Coefficients for the linear system (3.33) + (3.35) in the projection method.
double a_proj_[Nx + 2], b_proj_[Nx + 2];  // Only values from first_i_proj to (last_i_proj-1) are used
double d_proj_[Nx + 2];  // Values from first_i_proj to last_i_proj are used
double lambdax_proj_[Nx + 2];  // Only values from 0/1 to Nx are used
double d_lastRow = 100;

/*
 * Test Functions
 */

void writeCSV_d_proj();
void printLambdax();
void writeCSV_lambdax();
void writeCSV_ab_proj();
void printFirstRowU_test();

// To get the corresponding matrix, we choose to put (3.35) as the last row. Then we directly reduce the last row to make
// the matrix upper triangular, and we solve the system backwards.

// Our a and b are different from those in the article:
// a_this[i][j] = a_orig[i][j] - b_orig[i][j] / Dx,    b_this[i][j] = b_orig[i][j] / Dx.
// In our program here, a and b are chosen to be exactly the matrix entries in (3.33).
// Note: In the last version, we only store the inverse of a. We now store the values of a directly.

// d values represent the last row of the matrix, corresponding to (3.35). They store the original last matrix row and the
// intermediate values in the calculation of a and b. More importantly, d stores the linear transformation on c to directly
// make the matrix upper triangular.

// During the life of this program, the cache values of a, b, and d will remain the same once calculated.

// Refer to notes for details.

// Calculate a and b values
void fillCache_ab_proj() {
	for (int i = first_i_proj; i < last_i_proj; ++i) {  // Last index is last_i_proj - 1, see (3.33)
		double x = getCellCenterX(i), b = (*pB_fptr)(x) - pA;
		b_proj_[i] = b;
		a_proj_[i] = (*pBxDer_fptr)(x) * Dx - b;
	}
}

// Calculate d values. It is separated from the calculation of a and b for testing reasons
void fillCache_d_proj() {
	d_proj_[first_i_proj] = d_lastRow / a_proj_[first_i_proj];
	int last = last_i_proj - 1;
	// Loop over only up to (last_i_proj - 2), as each iteration calculates value of d[i+1]. d[last_i_proj] is defined alone as division is not needed.
	for (int i = first_i_proj; i < last; ++i)
		d_proj_[i + 1] = (d_lastRow - b_proj_[i] * d_proj_[i]) / a_proj_[i + 1];
	d_proj_[last_i_proj] = d_lastRow - b_proj_[last] * d_proj_[last];
	writeCSV_d_proj();
}

// Calculate the c_i's as in (3.34). Note that c_i changes in each time step.
// In this program, we do not have an array for c_i. Instead, we use lambda_x[] to hold the c_i values. Later after applying
// Gaussian elimination, lambda_x[] will naturally hold the lambda_x values as the solution of the linear system.
void calc_c_proj() {
	// Temporarily store int_pA^{pB(x_i)}\tilde{u}dp in c_i
	for (int i = first_i_proj; i <= last_i_proj; ++i) {  // Largest index is last_i_proj, b/c we need integral of u at x_{i+1}. See (3.33)
		double sum = 0;
		for (int j = 1; j <= Np; j++)
			sum += u_[i][j];
		sum *= getCellCenterDp(i);
		lambdax_proj_[i] = sum;
	}

	// Now calculate the real c_proj values
	for (int i = first_i_proj; i < last_i_proj; ++i)  // Only up to (last index - 1), as we need the (i+1)th element in the calculation of the ith element. See (3.34)
		lambdax_proj_[i] = lambdax_proj_[i + 1] - lambdax_proj_[i];
	lambdax_proj_[last_i_proj] = 0;  // This is in fact not necessary, since we directly assign values to c_i in the Gaussian elimination
	// writeCSV_lambdax();
}

// Perform Gaussian elimination to calculate lambda_x
void calcLambdax_proj() {
	calc_c_proj();
	double sum = 0;
	for (int i = first_i_proj; i < last_i_proj; ++i)  // Up to (last_i_proj-1)
		sum += lambdax_proj_[i] * d_proj_[i];
	lambdax_proj_[last_i_proj] = -sum / d_proj_[last_i_proj];
	for (int i = last_i_proj - 1; i >= first_i_proj; --i)
		lambdax_proj_[i] = (lambdax_proj_[i] - b_proj_[i] * lambdax_proj_[i + 1]) / a_proj_[i];
}

// Perform the projection method on uTilde to calculate u
void projU_orig() {
	calcLambdax_proj();
	for (int i = first_i_proj; i <= last_i_proj; ++i) {
		double lambdax_val = lambdax_proj_[i];
		for (int j = 1; j <= Np; ++j)
			u_[i][j] -= lambdax_val;
	}
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Test Functions
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void writeCSV_d_proj() {
	FILE *f = fopen("Output/d_proj.csv", "wb");
	for (int i = 1; i < Nx; ++i)
		fprintf(f, "%1.20e,", d_proj_[i]);
	fprintf(f, "%1.20e\n", d_proj_[Nx]);
	fclose(f);
}

void printLambdax() {
	printf("Lambdax values:  ");
	for (int i = 1; i < Nx; ++i)
		printf("%1.5e  ", lambdax_proj_[i]);
	printf("\n");
}

void writeCSV_lambdax() {
	FILE *f = fopen("Output/lambdax.csv", "wb");
	for (int i = 1; i < Nx; ++i)
		fprintf(f, "%1.20e,", lambdax_proj_[i]);
	fprintf(f, "%1.20e\n", lambdax_proj_[Nx]);
	fclose(f);
	ioWarning("writeCSV_lambdax()");
}

void writeCSV_ab_proj() {
	FILE *f = fopen("Output/ab_vec_proj.csv", "wb");
	int last = Nx - 1;
	for (int i = 1; i < last; ++i)
		fprintf(f, "%1.20e,", a_proj_[i]);
	fprintf(f, "%1.20e\n", a_proj_[last]);
	for (int i = 1; i < last; ++i)
		fprintf(f, "%1.20e,", b_proj_[i]);
	fprintf(f, "%1.20e\n", b_proj_[last]);
	fclose(f);
	printf("\n>> writeCSV_ab_proj() printed out something to CSV file!\n");
}

void printFirstRowU_test() {
	FILE *f = fopen("Output/u_last_row.csv", "wb");
	for (int i = 1; i < Nx; ++i)
		fprintf(f, "%1.20e,", u_[i][Np]);
	fprintf(f, "%1.20e", u_[Nx][Np]);
	printf("\n");
	fclose(f);
	printf(">> I/O: printFirstRowU_test printed something!\n");
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Set Projection Parameters and Fill Cache Values
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void setProjection() {
	fillCache_ab_proj();
	fillCache_d_proj();
	if (modelNo == 1)
		projU_fptr = &empty_fcn;
	else
		projU_fptr = &projU_orig;
}


#endif /* F_PROJECTION_H_ */
