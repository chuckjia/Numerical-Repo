/*
 * Mesh.h
 *
 *  Created on: Jun 20, 2017
 *      Author: chuckjia
 *
 *  This file contains methods needed for creating the mesh.
 */

#ifndef MESH_H_
#define MESH_H_
#include "Constants.h"

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Building the Mesh
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

/*
 * Initialize the space for the mesh
 */
double cellLeftX[numCellsXDir],  // Stores x-coords of left side of cell
cellRightX[numCellsXDir],  // Stores x-coords of right side of cell
cellTopP[numCellsXDir],  // Stores p-coords of top side of cell
cellBottP[numCellsXDir];  // Stores p-coords of bottom side of cell
double cellCenter[numCellsXDir][numCellsPDir][2];  // Stores coords of centers of all cells
double cellVol = Dx * Dp;

/*
 * Constructs the grid: calculate cellLeftX, cellRightX, cellTopP, cellBottP
 */
void buildGrid() {
	// Temporary values
	int tempInt;
	double tempDb;

	// x direction
	cellLeftX[0] = 0;
	cellRightX[0] = 0;
	tempInt = numCellsXDir - 1;
	for (int i = 1; i < tempInt; i++) {
		cellLeftX[i] = (i - 1) * Dx;
		cellRightX[i] = i * Dx;
	}
	tempDb = cellRightX[tempInt - 1];
	cellLeftX[tempInt] = tempDb;
	cellRightX[tempInt] = tempDb;

	// p direction
	cellTopP[0] = 0;
	cellBottP[0] = 0;
	tempInt = numCellsPDir - 1;
	for (int j = 1; j < tempInt; j++) {
		cellBottP[j] = (j - 1) * Dp;
		cellTopP[j] = j * Dp;
	}
	tempDb = cellTopP[tempInt - 1];
	cellBottP[tempInt] = tempDb;
	cellTopP[tempInt] = tempDb;
}

/*
 * Calculate the center coordinates for all cells and store them in the cellCenter
 */
void calcCellCenter() {
	int tempInt1 = numCellsXDir - 1, tempInt2 = numCellsPDir - 1;
	for (int i = 0; i < tempInt1; i++)
		for (int j = 0; j < tempInt2; j++) {
			cellCenter[i][j][0] = 0.5 * (cellLeftX[i] + cellRightX[i]);
			cellCenter[i][j][1] = 0.5 * (cellTopP[j] + cellBottP[j]);
		}
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Getters for cell coordinates and geometry
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

/*
 * Get the x coordinate on the left side of a cell
 */
double getCellLeftX(int i, int j) {
	return cellLeftX[i];
}

/*
 * Get the x coordinate on the right side of a cell
 */
double getCellRightX(int i, int j) {
	return cellRightX[i];
}

/*
 * Get the p coordinate on the top side of a cell
 */
double getCellTopP(int i, int j) {
	return cellTopP[j];
}

/*
 * Get the p coordinate on the bottom side of a cell
 */
double getCellBottP(int i, int j) {
	return cellBottP[j];
}

/*
 * Get the x coordinate of the center of a cell
 */
double getCenterX(int i, int j) {
	return cellCenter[i][j][0];
}

/*
 * Get the p coordinate of the center of a cell
 */
double getCenterP(int i, int j) {
	return cellCenter[i][j][1];
}

/*
 * Get the volume of a cell
 */
double getVol(int i, int j) {
	if (i == 0 || j == 0 || i == numGridPtsXDir || j == numGridPtsPDir)
		return 0;
	return cellVol;
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Wrapper function to build the Mesh
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void buildMesh() {
	buildGrid();
	calcCellCenter();
}

#endif /* MESH_H_ */
