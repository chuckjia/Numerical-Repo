/*
 * TestFcns.h
 *
 *  Created on: Jul 17, 2017
 *      Author: chuckjia
 */

#ifndef TESTFCNS_H_
#define TESTFCNS_H_
#include "FV.h"

void printMat(int numRows, int numCols, double mat[numRows][numCols]) {
	for (int i = 0; i < numRows; i++) {
		for (int j = 0; j < numCols; j++) {
			printf("%1.2f ", mat[i][j]);
		}
		printf("\n");
	}
}

void printFluxGG() {
	printf("\nGG = \n");
	for (int i = 1; i < lastGhostIndexX; i++) {
		for (int j = 0; j < lastGhostIndexP; j++) {
			printf("%1.7f ", GG[i][j]);
		}
		printf("\n");
	}
}

void printFluxFF() {
	printf("\nFF = \n");
	for (int i = 0; i < lastGhostIndexX; i++) {
		for (int j = 1; j < lastGhostIndexP; j++) {
			printf("%1.7f ", FF[i][j]);
		}
		printf("\n");
	}
}

#endif /* TESTFCNS_H_ */
