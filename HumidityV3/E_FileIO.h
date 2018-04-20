/*
 * FileIO.h
 *
 *  Created on: Apr 6, 2018
 *      Author: chuckjia
 */

#ifndef E_FILEIO_H_
#define E_FILEIO_H_

#include "D_Mesh.h"

void printMeshGridToFile(const char* filename_x, const char* filename_p) {
	FILE *fx = fopen(filename_x, "wb");
	FILE *fp = fopen(filename_p, "wb");

	int last_j = numGridPtsP - 1;

	for (int i = 0; i < numGridPtsX; ++i) {
		double gridX = meshGridX[i];
		for (int j = 0; j < last_j; ++j) {
			double gridP = meshGridP[i][j];
			fprintf(fx, "%1.20e,", gridX);
			fprintf(fp, "%1.20e,", gridP);
		}
		fprintf(fx, "%1.20e\n", gridX);
		fprintf(fp, "%1.20e\n", meshGridP[i][last_j]);
	}

	fclose(fx); fclose(fp);
}

void printMeshGridToFile() {
	printMeshGridToFile("Output/MeshGrid_X.csv", "Output/MeshGrid_P.csv");
}

void printMatrixToFile(double mat[numCellsX][numCellsP], const char* filename) {
	FILE *f = fopen(filename, "wb");
	int last_j = numCellsP - 1;
	for (int i = 0; i < numCellsX; ++i) {
		for (int j = 0; j < last_j; ++j)
			fprintf(f, "%1.20e,", mat[i][j]);
		fprintf(f, "%1.20e\n", mat[i][last_j]);
	}
	fclose(f);
}

void printMatrixToFile(double mat[numCellsX][numCellsP], const char* prefix, int step) {
	char filename[30];
	sprintf(filename, "%s_%d.csv", prefix, step);
	printMatrixToFile(mat, filename);
}

void printSolnToFile(int tt) {
	char filename[30];
	sprintf(filename, "MovieFrames/T_%d.csv", tt);
	printMatrixToFile(T_sl, filename);

	sprintf(filename, "MovieFrames/q_%d.csv", tt);
	printMatrixToFile(q_sl, filename);

	sprintf(filename, "MovieFrames/u_%d.csv", tt);
	printMatrixToFile(u_sl, filename);

	sprintf(filename, "MovieFrames/w_%d.csv", tt);
	printMatrixToFile(w_sl, filename);
}

// Print cell center coordinates to file
void printCellCentersToFile() {
	// Print x-coordinates
	FILE *f = fopen("Output/CellCenters_X.csv", "wb");
	int last_j = numCellsP - 1;
	for (int i = 0; i < numCellsX; ++i) {
		double centerX = cellCenterX[i];
		for (int j = 0; j < last_j; ++j)
			fprintf(f, "%1.20e,", centerX);
		fprintf(f, "%1.20e\n", centerX);
	}
	fclose(f);

	// Print p-coordinates
	printMatrixToFile(cellCenterP, "Output/CellCenters_P.csv");
}



#endif /* E_FILEIO_H_ */
