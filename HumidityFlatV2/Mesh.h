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
double cellSidesX[numCellsX][2];
// Stores p-coords of the bottom (0) and top (1) sides of a cell
double cellSidesP[numCellsP][2];
double cellCenters[numCellsX][numCellsP][2];  // Stores coords of centers of all cells
double cellVol = Dx * Dp; // The volume of all cells (in this model, all have same volume)
double cellVolInv = DxInv * DpInv;

/*
 * Constructs the grid: calculate cellSidesX and cellSidesP
 */
void buildGrid() {
	// x direction
	cellSidesX[0][0] = x0;
	cellSidesX[0][1] = x0;
	for (int i = 1; i < lastIndexX; i++) {
		cellSidesX[i][0] = x0 + (i - 1) * Dx;
		cellSidesX[i][1] = x0 + i * Dx;
	}
	cellSidesX[lastIndexX][0] = xL;
	cellSidesX[lastIndexX][1] = xL;

	// p direction
	cellSidesP[0][0] = pA;
	cellSidesP[0][1] = pA;
	for (int j = 1; j < lastIndexP; j++) {
		cellSidesP[j][0] = pA + (j - 1) * Dp;
		cellSidesP[j][1] = pA + j * Dp;
	}
	cellSidesP[lastIndexP][0] = pB;
	cellSidesP[lastIndexP][1] = pB;
}

/*
 * Calculate the center coordinates for all cells and store them in cellCenters
 */
void calcCellCenter() {
	for (int i = 0; i < numCellsX; i++)
		for (int j = 0; j < numCellsP; j++) {
			cellCenters[i][j][0] = 0.5 * (cellSidesX[i][0] + cellSidesX[i][1]);
			cellCenters[i][j][1] = 0.5 * (cellSidesP[j][0] + cellSidesP[j][1]);
		}
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Getters for cell coordinates and geometry
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

/*
 * Get the x coordinate on the left side of a cell
 */
double getCellLeftX(int i, int j) {
	return cellSidesX[i][0];
}

/*
 * Get the x coordinate on the right side of a cell
 */
double getCellRightX(int i, int j) {
	return cellSidesX[i][1];
}

/*
 * Get the p coordinate on the bottom side of a cell
 */
double getCellBottP(int i, int j) {
	return cellSidesP[j][0];
}

/*
 * Get the p coordinate on the top side of a cell
 */
double getCellTopP(int i, int j) {
	return cellSidesP[j][1];
}

/*
 * Get the x coordinate of the center of a cell
 */
double getCellCenterX(int i, int j) {
	return cellCenters[i][j][0];
}

/*
 * Get the p coordinate of the center of a cell
 */
double getCellCenterP(int i, int j) {
	return cellCenters[i][j][1];
}

/*
 * Get the coordinate of the center of a cell and pass to a pointer
 */
/*void getCellCenter(double ans[2], int i, int j) {
	ans[0] = cellCenters[i][j][0];
	ans[1] = cellCenters[i][j][1];
}*/

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
	// Measure execution time
	clock_t startTime, endTime;
	startTime = clock();

	// Main body
	buildGrid();
	calcCellCenter();

	// Measure execution time
	endTime = clock();
	double cpuTimeUsed = (double) (endTime - startTime) * 0.001;
	printf("\n- Mesh successfully built. (%1.2fms)\n", cpuTimeUsed);
}

#endif /* MESH_H_ */
