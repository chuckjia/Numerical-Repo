/*
 * E_FileIO.h
 *
 *  Created on: Apr 6, 2018
 *      Author: chuckjia
 *
 *  This file includes functions that perform file I/O operations.
 */

#ifndef E_FILEIO_H_
#define E_FILEIO_H_

#include "D_Mesh.h"

// Print the mesh grid points to file
// The x- and p-coordinates of the grid points are stored in two different files, with names filenameX and filenameP.
void printMeshGridToFile(const char* filenameX, const char* filenameP) {
	FILE *fx = fopen(filenameX, "wb");
	FILE *fp = fopen(filenameP, "wb");

	int last_j = numGridPtP - 1;  // Used to print the last entry of each row, which does not have commas following

	for (int i = 0; i < numGridPtX; ++i) {
		double gridX = meshGridX_[i];
		for (int j = 0; j < last_j; ++j) {
			double gridP = meshGridP_[i][j];
			fprintf(fx, "%1.20e,", gridX);
			fprintf(fp, "%1.20e,", gridP);
		}
		fprintf(fx, "%1.20e\n", gridX);
		fprintf(fp, "%1.20e\n", meshGridP_[i][last_j]);
	}

	fclose(fx); fclose(fp);
}

// Print the mesh grid points to files with standard names
void printMeshGridToFile() {
	printMeshGridToFile("Output/MeshGrid_X.csv", "Output/MeshGrid_P.csv");
}

// Print matrices with the same size as the solution matrices to file
void printMatrixToFile(double mat[numCellX][numCellP], const char* filename) {
	FILE *f = fopen(filename, "wb");
	int last_j = numCellP - 1;
	for (int i = 0; i < numCellX; ++i) {
		for (int j = 0; j < last_j; ++j)
			fprintf(f, "%1.20e,", mat[i][j]);
		fprintf(f, "%1.20e\n", mat[i][last_j]);
	}
	fclose(f);
}

// Print solutions to file. File names are appended with time step numbers
void printSolnToFile(int tt) {
	char filename[30];
	sprintf(filename, "MovieFrames/T_%d.csv", tt);
	printMatrixToFile(T_, filename);

	sprintf(filename, "MovieFrames/q_%d.csv", tt);
	printMatrixToFile(q_, filename);

	sprintf(filename, "MovieFrames/u_%d.csv", tt);
	printMatrixToFile(u_, filename);

	sprintf(filename, "MovieFrames/w_%d.csv", tt);
	printMatrixToFile(w_, filename);
}

// Print cell center coordinates to file
void printCellCentersToFile() {
	// Print x-coordinates
	FILE *f = fopen("Output/CellCenters_X.csv", "wb");
	int last_j = numCellP - 1;
	for (int i = 0; i < numCellX; ++i) {
		double centerX = cellCenterX_[i];
		for (int j = 0; j < last_j; ++j)
			fprintf(f, "%1.20e,", centerX);
		fprintf(f, "%1.20e\n", centerX);
	}
	fclose(f);

	// Print p-coordinates
	printMatrixToFile(cellCenterP_, "Output/CellCenters_P.csv");
}



#endif /* E_FILEIO_H_ */
