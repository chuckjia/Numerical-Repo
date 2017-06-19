/*
 * TestFunc.h
 *
 *  Created on: Jun 18, 2017
 *      Author: chuckjia
 */

#ifndef TESTFCNS_H_
#define TESTFCNS_H_

void printGrid(int numRows, int numCols, double mat[numRows][numCols]) {
	for (int i = 0; i < numRows; i++) {
		for (int j = 0; j < numCols; j++) {
			printf("%5.2f ", mat[i][j]);
		}
		printf("\n");
	}
}

#endif /* TESTFCNS_H_ */
