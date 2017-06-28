/*
 * Mesh.h
 *
 *  Created on: Jun 10, 2017
 *      Author: chuckjia
 */

#ifndef MESH_H_
#define MESH_H_

#include "Constants.h"

double DpArr[numGridPtsXDir];  // DpArr[i] stores the Delta p value at x[i + 1/2][j]

/*
 * meshGridY[i][j] represents the y coordinate of the grid
 * point (i + 1/2, j + 1/2)
 * meshGrid is global
 */
double meshGridP[numGridPtsXDir][numGridPtsPDir];

// The following are arrays storing the geometry of cells
double cellVols[numCellsX][numCellsP];  // Stores the volumes of all cells
double cellCenters[numCellsX][numCellsP][2];  // Stores the centers of all cells

/*
 * Calculate and store all the Delta p values. DpArr[i] represents the Delta p
 * value at the right side of the cell (i, j). The right most flat control volume
 * is not included.
 */

void calcDp() {
	for (int i = 0; i < numGridPtsXDir; i++)
		DpArr[i] = Dp(x0 + i * Dx);
}

/*
 * Getter function for the Delta p values
 */

double getDp(int i) {
	if (i == numGridPtsXDir)
		return DpArr[i - 1];
	return DpArr[i];
}

/*
 * buildGrid constructs the grid meshGridY
 */

void buildGrid() {
	for (int i = 0; i < numGridPtsXDir; i++) {
		double DpVal = DpArr[i];
		for (int j = 0; j < numGridPtsPDir; j++)
			meshGridP[i][j] = j * DpVal;
	}
}

/*
 * Following functions are getters for coordinates of each cell
 */

double getCellRightX(int i, int j) {
	if (i == numGridPtsXDir)
		return xL;
	return i * Dx;
}

double getCellLeftX(int i, int j) {
	if (i == 0)
		return x0;
	return (i - 1) * Dx;
}

double getTopRightP(int i, int j) {
	if (i == numGridPtsXDir)
		i--;
	if (j == numGridPtsPDir)
		j--;
	return meshGridP[i][j];
}

double getTopLeftP(int i, int j) {
	if (i == 0)
		i++;
	if (j == numGridPtsPDir)
		j--;
	return meshGridP[i - 1][j];
}

double getBottLeftP(int i, int j) {
	if (i == 0)
		i++;
	if (j == 0)
		j++;
	return meshGridP[i - 1][j - 1];
}

double getBottRightP(int i, int j) {
	if (i == numGridPtsXDir)
		i--;
	if (j == 0)
		j++;
	return meshGridP[i][j - 1];
}

/*
 * This function calculates the volume of a cell
 */

double calcVol(double topLeft[2], double topRight[2], double bottLeft[2],
		double bottRight[2]) {
	double h = topRight[0] - topLeft[0];
	return 0.5 * h * (topLeft[1] - bottLeft[1] + topRight[1] - bottRight[1]);
}

/*
 * This function calculates the bary-center of a given cell
 */

void calcBaryCenter(double center[2], double topLeft[2], double topRight[2],
		double bottLeft[2], double bottRight[2]) {
	// South-East triangle
	double xbar1 = (bottLeft[0] + bottRight[0] + topRight[0]) / 3.0;
	double pbar1 = (bottLeft[1] + bottRight[1] + topRight[1]) / 3.0;
	// North-East triangle
	double xbar3 = (topLeft[0] + bottRight[0] + topRight[0]) / 3.0;
	double pbar3 = (topLeft[1] + bottRight[1] + topRight[1]) / 3.0;
	// North-West triangle
	double xbar2 = (topLeft[0] + bottLeft[0] + topRight[0]) / 3.0;
	double pbar2 = (topLeft[1] + bottLeft[1] + topRight[1]) / 3.0;
	// South-West triangle
	double xbar4 = (topLeft[0] + bottLeft[0] + bottRight[0]) / 3.0;
	double pbar4 = (topLeft[1] + bottLeft[1] + bottRight[1]) / 3.0;
	// Calculate the intersection of two diagonal lines
	double denom = (xbar1 - xbar2) * (pbar3 - pbar4) - (pbar1 - pbar2) * (xbar3 - xbar4);
	center[0] = ((xbar1 * pbar2 - pbar1 * xbar2) * (xbar3 - xbar4) - (xbar1 - xbar2) * (
			xbar3 * pbar4 - pbar3 * xbar4)) / denom;
	center[1] = ((xbar1 * pbar2 - pbar1 * xbar2) * (pbar3 - pbar4) - (pbar1 - pbar2) * (
			xbar3 * pbar4 - pbar3 * xbar4)) / denom;
}

/*
 * Calculate the volume and the bary-center for each cell
 */

void calcGeometry() {
	// Calculate for all normal (inner) cells
	int rightmostXInd = numCellsX - 1, upmostPInd = numCellsP - 1;
	for (int i = 1; i < rightmostXInd; i++)
		for (int j = 1; j < upmostPInd; j++) {
			// Cell corner coordinates
			double xLeft = getCellLeftX(i, j), xRight = getCellRightX(i, j);
			double topLeft[2] = {xLeft, getTopLeftP(i, j)};
			double topRight[2] = {xRight, getTopRightP(i, j)};
			double bottLeft[2] = {xLeft, getBottLeftP(i, j)};
			double bottRight[2] = {xRight, getBottRightP(i, j)};
			// Cell volume
			cellVols[i][j] = calcVol(topLeft, topRight, bottLeft, bottRight);
			// Cell bary-center
			double center[2];
			calcBaryCenter(center, topLeft, topRight, bottLeft, bottRight);
			cellCenters[i][j][0] = center[0];
			cellCenters[i][j][1] = center[1];
		}
	// The following are for flat control volumes
	// When i = 0 and i = numCellsXDir - 1 (rightmostX)
	int iLeft = 0, iRight = rightmostXInd;
	for (int j = 1; j < upmostPInd; j++) {
		// Cell volume
		cellVols[iLeft][j] = 0;
		cellVols[iRight][j] = 0;
		// Cell center
		cellCenters[iLeft][j][0] = x0;
		cellCenters[iLeft][j][1] = 0.5 * (getTopRightP(iLeft, j) + getBottRightP(iLeft, j));
		cellCenters[iRight][j][0] = xL;
		cellCenters[iRight][j][1] = 0.5 * (getTopRightP(iRight, j) + getBottRightP(iRight, j));
	}
	// When j = 0 and j = numCellsPDir - 1 (upmostP)
	int jBott = 0, jTop = upmostPInd;
	for (int i = 1; i < rightmostXInd; i++) {
		// Cell volume
		cellVols[i][jBott] = 0;
		cellVols[i][jTop] = 0;
		// Cell center
		double xCenter = 0.5 * (getCellLeftX(i, jBott) + getCellRightX(i, jBott));
		cellCenters[i][jBott][0] = xCenter;
		cellCenters[i][jBott][1] = pA;
		cellCenters[i][jTop][0] = xCenter;
		cellCenters[i][jTop][1] = 0.5 * (getBottLeftP(i, jTop) + getBottRightP(i, jTop));
	}
	// Corners
	cellVols[0][0] = 0;
	cellVols[0][upmostPInd] = 0;
	cellVols[rightmostXInd][0] = 0;
	cellVols[rightmostXInd][upmostPInd] = 0;
	cellCenters[0][0][0] = x0;
	cellCenters[0][0][1] = pA;
	cellCenters[0][upmostPInd][0] = x0;
	cellCenters[0][upmostPInd][1] = pB(x0);
	cellCenters[rightmostXInd][0][0] = xL;
	cellCenters[rightmostXInd][0][1] = pA;
	cellCenters[rightmostXInd][upmostPInd][0] = xL;
	cellCenters[rightmostXInd][upmostPInd][1] = pB(xL);
}

double getVol(int i, int j) {
	return cellVols[i][j];
}

double getCellCenterX(int i, int j) {
	return cellCenters[i][j][0];
}

double getCellCenterP(int i, int j) {
	return cellCenters[i][j][1];
}

void buildMesh() {
	calcDp();
	buildGrid();
	calcGeometry();
}

#endif /* MESH_H_ */
