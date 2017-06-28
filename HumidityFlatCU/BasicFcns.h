/*
 * BasicFcns.h
 *
 *  Created on: Jun 10, 2017
 *      Author: chuckjia
 *
 *  This file contains mathematical functions that are independent from
 *  model settings.
 */

#ifndef BASICFCNS_H_
#define BASICFCNS_H_
#include <stdio.h>
#include <math.h>
#include <time.h>

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Basic Math Functions
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

/*
 * The method realizes the mathematical function min(x, y)
 */
double min(double x, double y) {
	if (x < y)
		return x;
	return y;
}

/*
 * The method realizes the mathematical function max(x, y)
 */
double max(double x, double y) {
	if (x > y)
		return x;
	return y;
}

/*
 * The method realizes the mathematical function min(x, y, z)
 */
double min3(double x, double y, double z) {
	if (x < y) {
		if (x < z)
			return x;
		return z;
	}
	// When y <= x
	if (y < z)
		return y;
	return z;
}

/*
 * The method realizes the mathematical function max(x, y, z)
 */
double max3(double x, double y, double z) {
	if (x > y) {
		if (x > z)
			return x;
		return z;
	}
	// When y >= x
	if (y > z)
		return y;
	return z;
}

/*
 * The method realizes the mathematical function minmod(x, y, z)
 */
double minmod3(double x, double y, double z) {
	if (x > 0 && y > 0 && z > 0)
		return min3(x, y, z);
	if (x < 0 && y < 0 && z < 0)
		return max3(x, y, z);
	return 0;
}

/*
 * The method realizes the sign function sgn(x, y, z)
 */
double sign(double x) {
	if (x > 0)
		return 1;
	if (x < 0)
		return -1;
	return 0;
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Methods For Testing Purposes
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

/*
 * Print out a 2D array/matrix of double-precision entries
 */
void printMatrix(int numRows, int numCols, double mat[numRows][numCols]) {
	for (int i = 0; i < numRows; i++) {
		for (int j = 0; j < numCols; j++) {
			printf("%5.2f ", mat[i][j]);
		}
		printf("\n");
	}
}

/*
 * Print out a 3D array/matrix of double-precision entries with the last
 * dimension having fixed length 2. That is, the matrix is essentially a 2D
 * array of pairs of double-precision numbers.
 */
void printMatrix2DTimes2(int numRows, int numCols, double mat[numRows][numCols][2]) {
	for (int i = 0; i < numRows; i++) {
		for (int j = 0; j < numCols; j++) {
			printf("(%0.2f, %0.2f) ", mat[i][j][0], mat[i][j][1]);
		}
		printf("\n");
	}
}

/*
 * This function works the same with printMatrix3D2(), except that it prints out
 * only the first component of each entry of the 2D array of pairs.
 */
void printMatrix2DTimes2_Comp1(int numRows, int numCols, double mat[numRows][numCols][2]) {
	for (int i = 0; i < numRows; i++) {
		for (int j = 0; j < numCols; j++) {
			double val = mat[i][j][0];
			if (val >= 0)
				printf("%3.3f ", val);
			else
				printf("%3.3f ", -val);
		}
		printf("\n");
	}
}

#endif /* BASICFCNS_H_ */
