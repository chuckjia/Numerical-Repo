/*
 * TestFcns.h
 *
 *  Created on: Jun 20, 2017
 *      Author: chuckjia
 *
 *  This header file contains functions created for testing purposes.
 */

#ifndef TESTFCNS_H_
#define TESTFCNS_H_
#include <stdio.h>
#include <math.h>

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
void printMatrix3D2(int numRows, int numCols, double mat[numRows][numCols][2]) {
	for (int i = 0; i < numRows; i++) {
		for (int j = 0; j < numCols; j++) {
			printf("(%5.1f, %5.1f) ", mat[i][j][0], mat[i][j][1]);
		}
		printf("\n");
	}
}

/*
 * This function works the same with printMatrix3D2(), except that it prints out
 * only the first component of each entry of the 2D array of pairs.
 */
void printMatrix3D_Comp1(int numRows, int numCols, double mat[numRows][numCols][2]) {
	for (int i = 0; i < numRows; i++) {
		for (int j = 0; j < numCols; j++) {
			printf("%5.1f ", mat[i][j][0]);
		}
		printf("\n");
	}
}

#endif /* TESTFCNS_H_ */
