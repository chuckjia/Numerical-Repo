/*
 * Testing.h
 *
 *  Created on: Oct 16, 2017
 *      Author: chuckjia
 */

#ifndef TESTING_H_
#define TESTING_H_
#include "Projection.h"

void printMeshToFile() {
	FILE *fGridX = fopen("gridX.txt", "wb");
	FILE *fGridP = fopen("gridP.txt", "wb");
	for (int i = 0; i < numGridPtsX; ++i)
		for (int j = 0; j < numGridPtsP; ++j) {
			fprintf(fGridX, "%f ", getCellLeftX(i, j));
			fprintf(fGridP, "%f ", getCellBottLeftP(i, j));
		}
	fclose(fGridX);
	fclose(fGridP);

	FILE *fCellCenterX = fopen("cellCenterX.txt", "wb");
	FILE *fCellCenterP = fopen("cellCenterP.txt", "wb");
	FILE *fCellVol = fopen("cellVol.txt", "wb");
	for (int i = 0; i < numCellsX; ++i)
		for (int j = 0; j < numCellsP; ++j) {
			fprintf(fCellCenterX, "%f ", getCellCenterX(i, j));
			fprintf(fCellCenterP, "%f ", getCellCenterP(i, j));
			fprintf(fCellVol, "%f ", getCellVol(i, j));
		}
	fclose(fCellCenterX);
	fclose(fCellCenterP);
	fclose(fCellVol);
}

void printPar() {
	FILE *f = fopen("par.txt", "wb");
	fprintf(f, "%f %f %f %f %f %d %d",
			x0, xf, pA, (*pB_fcnPtr)(x0), (*pB_fcnPtr)(xf), Nx, Np);
	fclose(f);
}

void printDiagnostics() {
	printf(">> DIAGNOSTIC INFORMATION\n");
	printf("\nParameters: \n");

	printf("\n - Domain geometry\n");
	printf("    [1] x0 = %1.2f, xf = %1.2f \n", x0, xf);
	printf("    [2] pA = %1.2f, pB(x0) = %1.2f, pB[(x0+xf)/2] = %1.2f, pB(xf) = %1.2f\n",
			pA, (*pB_fcnPtr)(x0), (*pB_fcnPtr)(0.5 * (x0 + xf)), (*pB_fcnPtr)(xf));

	printf("\n - Mesh specifications\n");
	printf("    [1] Nx = %d, Np = %d\n", Nx, Np);
	printf("    [2] Dx = %1.2f\n", Dx);
	printf("    [3] Dp(x0) = %1.2f, Dp[(x0+xf)/2] = %1.2f, Dp(xf) = %1.2f",
			getCellLeftDp(1), getCellLeftDp(lastRealIndexX / 2), getCellLeftDp(lastRealIndexX));
}

#endif /* TESTING_H_ */
