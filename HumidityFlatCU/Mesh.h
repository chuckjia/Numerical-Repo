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
// Stores x-coords of the left (0) and right (1) sides of a cell
double cellSideX[numCellsX][2];
// Stores p-coords of the bottom (0) and top (1) sides of a cell
double cellSideP[numCellsP][2];
double cellCenter[numCellsX][numCellsP][2];  // Stores coords of centers of all cells
double cellVol = Dx * Dp;

/*
 * Constructs the grid: calculate cellLeftX, cellRightX, cellTopP, cellBottP
 */
void buildGrid() {
	// x direction
	cellSideX[0][0] = x0;
	cellSideX[0][1] = x0;
	for (int i = 1; i < lastIndexX; i++) {
		cellSideX[i][0] = x0 + (i - 1) * Dx;
		cellSideX[i][1] = x0 + i * Dx;
	}
	cellSideX[lastIndexX][0] = xL;
	cellSideX[lastIndexX][1] = xL;

	// p direction
	cellSideP[0][0] = 0;
	cellSideP[0][1] = 0;
	for (int j = 1; j < lastIndexP; j++) {
		cellSideP[j][0] = pA + (j - 1) * Dp;
		cellSideP[j][1] = pA + j * Dp;
	}
	cellSideP[lastIndexP][0] = pB;
	cellSideP[lastIndexP][1] = pB;
}

/*
 * Calculate the center coordinates for all cells and store them in the cellCenter
 */
void calcCellCenter() {
	for (int i = 0; i < numCellsX; i++)
		for (int j = 0; j < numCellsP; j++) {
			cellCenter[i][j][0] = 0.5 * (cellSideX[i][0] + cellSideX[i][1]);
			cellCenter[i][j][1] = 0.5 * (cellSideP[j][0] + cellSideP[j][1]);
		}
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Getters for cell coordinates and geometry
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

/*
 * Get the x coordinate on the left side of a cell
 */
double getCellLeftX(int i, int j) {
	return cellSideX[i][0];
}

/*
 * Get the x coordinate on the right side of a cell
 */
double getCellRightX(int i, int j) {
	return cellSideX[i][1];
}

/*
 * Get the p coordinate on the bottom side of a cell
 */
double getCellBottP(int i, int j) {
	return cellSideP[j][0];
}

/*
 * Get the p coordinate on the top side of a cell
 */
double getCellTopP(int i, int j) {
	return cellSideP[j][1];
}

/*
 * Get the x coordinate of the center of a cell
 */
double getCellCenterX(int i, int j) {
	return cellCenter[i][j][0];
}

/*
 * Get the p coordinate of the center of a cell
 */
double getCellCenterP(int i, int j) {
	return cellCenter[i][j][1];
}

/*
 * Get the coordinate of the center of a cell and pass to a pointer
 */
void getCenter(double ans[2], int i, int j) {
	ans[0] = cellCenter[i][j][0];
	ans[1] = cellCenter[i][j][1];
}

/*
 * Get the volume of a cell
 */
double getVol(int i, int j) {
	if (i == 0 || j == 0 || i == lastIndexX || j == lastIndexP)
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
