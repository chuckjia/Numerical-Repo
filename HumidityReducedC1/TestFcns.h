/*
 * TestFcns.h
 *
 *  Created on: Jun 10, 2017
 *      Author: chuckjia
 */

#ifndef TESTFCNS_H_
#define TESTFCNS_H_
#include <stdio.h>
#include <math.h>

void printMatrix(int numRows, int numCols, double mat[numRows][numCols]) {
	for (int i = 0; i < numRows; i++) {
		for (int j = 0; j < numCols; j++) {
			printf("%5.2f ", mat[i][j]);
		}
		printf("\n");
	}
}

void printMatrix3D(int numRows, int numCols, double mat[numRows][numCols][2]) {
	for (int i = 0; i < numRows; i++) {
		for (int j = 0; j < numCols; j++) {
			printf("(%5.1f, %5.1f) ", mat[i][j][0], mat[i][j][1]);
		}
		printf("\n");
	}
}

void printMatrix3D_1(int numRows, int numCols, double mat[numRows][numCols][2]) {
	for (int i = 0; i < numRows; i++) {
		for (int j = 0; j < numCols; j++) {
			printf("%5.1f ", mat[i][j][0]);
		}
		printf("\n");
	}
}

#endif /* TESTFCNS_H_ */
