/*
 * Mesh.h
 *
 *  Created on: Oct 11, 2017
 *      Author: chuckjia
 */

#ifndef MESH_H_
#define MESH_H_
#include "Models.h"

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Mesh Constants: Auto-generated
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

const int Nx = numDivisions;  // Number of space divisions in the x direction
const int Np = numDivisions;  // Number of space divisions in the p direction

// Following constants are set up for code readability
double NxInv = 1. / Nx, NpInv = 1. / Np;
const int numCellsX = Nx + 2;  // Number of cells in the x direction
const int numCellsP = Np + 2;  // Number of cells in the p direction
int lastRealIndexX = Nx;  // Largest x-index of all non-ghost cells
int lastRealIndexP = Np;  // Largest p-index of all non-ghost cells
int lastGhostIndexX = Nx + 1;  // x-index of ghost cells on the right side of domain
int lastGhostIndexP = Np + 1;  // p-index of ghost cells on the top side of domain
const int numGridPtsX = numCellsX + 1;  // Num of grid points in x direction (including ghosts)
const int numGridPtsP = numCellsP + 1;  // Num of grid points in p direction (including ghosts)
const int NxPlusOne = Nx + 1;
const int NpPlusOne = Np + 1;

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Grid Points
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

double Dx, DxInv;  // x step size, defined in calcGridPts()
double cellLeftDpVec[numGridPtsX];  // Stores the Dp values on cell left sides

// The following two arrays store the grid point coordinates. meshgridX stores x coord on cell
// left sides, and meshgridP stores p coord on cell bottom-left sides
double meshgridX[numGridPtsX];  // Stores x_{i-1/2,j}
double meshgridP[numGridPtsX][numGridPtsP];  // Stores p_{i-1/2,j-1/2}

// Calculate and store grid points coordinates in cache
void calcGridPts() {
	Dx = (xf - x0) / Nx;
	DxInv = Nx / (xf - x0);
	for (int i = 0; i < numGridPtsX; ++i) {
		double x = x0 + (i - 1) * Dx;
		meshgridX[i] = x;
		double Dp = ((*pB_fcnPtr)(x) - pA) * NpInv;
		cellLeftDpVec[i] = Dp;
		for (int j = 0; j < numGridPtsP; ++j)
			meshgridP[i][j] = pA + (j - 1) * Dp;
	}
}

// Getters for the Dp values

double getCellLeftDp(int i, int j) {
	return cellLeftDpVec[i];
}

double getCellLeftDp(int i) {
	return cellLeftDpVec[i];
}

double getCellRightDp(int i, int j) {
	return cellLeftDpVec[i + 1];
}

double getCellRightDp(int i) {
	return cellLeftDpVec[i + 1];
}

// Getters for x coordinates on grid points

double getCellLeftX(int i, int j) {
	return meshgridX[i];
}

double getCellLeftX(int i) {
	return meshgridX[i];
}

double getCellRightX(int i, int j) {
	return meshgridX[i + 1];
}

double getCellRightX(int i) {
	return meshgridX[i + 1];
}

// Getters for p coordinates on grid points

double getCellBottLeftP(int i, int j) {
	return meshgridP[i][j];
}

double getCellBottRightP(int i, int j) {
	return meshgridP[i + 1][j];
}

double getCellTopLeftP(int i, int j) {
	return meshgridP[i][j + 1];
}

double getCellTopRightP(int i, int j) {
	return meshgridP[i + 1][j + 1];
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Cell Centers
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Store the coordinates for cell barycenters
double cellCenterX[numCellsX], cellCenterP[numCellsX][numCellsP];

// Calculate cell barycenters. Require cells to be trapezoids
void calcBaryCenters() {
	double oneThird = 1. / 3;
	double h = Dx / 3.;  // h = (x2_Outer - x1_Outer) / 3
	for (int i = 0; i < numCellsX; ++i) {
		double x1_Outer = getCellLeftX(i), x2_Outer = getCellRightX(i);
		double x1_Inner = (2 * x1_Outer + x2_Outer) * oneThird;
		double a = getCellRightDp(i) * oneThird,  // a = (p3_Outer - p2_Outer) / 3
				b = getCellLeftDp(i) * oneThird,  // b = (p4_Outer - p1_Outer) / 3
				inv_ab_sum = 1. / (a + b);
		cellCenterX[i] = x1_Inner + h * a * inv_ab_sum;
		for (int j = 0; j < numCellsP; ++j) {
			double p1_Outer = getCellBottLeftP(i, j), p2_Outer = getCellBottRightP(i, j),
					p3_Outer = getCellTopRightP(i, j), p4_Outer = getCellTopLeftP(i, j);
			double p1_Inner = (p1_Outer + p2_Outer + p4_Outer) * oneThird;
			double bLeft = (p4_Outer - p3_Outer) * oneThird;
			cellCenterP[i][j] = p1_Inner + (b - bLeft) * a * inv_ab_sum;
		}
	}
}

// Getters for coordinates of cell centers

double getCellCenterX(int i, int j) {
	return cellCenterX[i];
}

double getCellCenterX(int i) {
	return cellCenterX[i];
}

double getCellCenterP(int i, int j) {
	return cellCenterP[i][j];
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Cell Side Lengths
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Stores the length of cell bottom sides
double cellBottSideLenVec[numCellsX][numCellsP];

// Calculate all cell sides length: calculate bottom side length of cell (i, j)
void calcCellSideLen() {
	for (int i = 0; i < numCellsX; ++i) {
		double xDiffSq = pow(getCellRightX(i) - getCellLeftX(i), 2);
		for (int j = 0; j < numCellsP; ++j) {
			cellBottSideLenVec[i][j] = sqrt(xDiffSq +
					pow(getCellBottRightP(i, j) - getCellBottLeftP(i, j), 2));
		}
	}
}

// Getters for cell side lengths

double getCellBottSideLen(int i, int j) {
	return cellBottSideLenVec[i][j];
}

double getCellTopSideLen(int i, int j) {
	return cellBottSideLenVec[i][j + 1];
}

double getCellLeftSideLen(int i, int j) {
	return cellLeftDpVec[i];
}

double getCellLeftSideLen(int i) {
	return cellLeftDpVec[i];
}

double getCellRightSideLen(int i, int j) {
	return cellLeftDpVec[i + 1];
}

double getCellRightSideLen(int i) {
	return cellLeftDpVec[i + 1];
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Cell Volumes
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Store the volumes of all cells. Made 2D to accommodate flat control volumes
double cellVol[numCellsX][numCellsP];

// Calculate cell volumes. Require cells to be trapezoids
void calcCellVol() {
	for (int i = 0; i < numCellsX; ++i)
		for (int j = 0; j < numCellsP; ++j)
			cellVol[i][j] = 0.5 * (getCellLeftDp(i) + getCellRightDp(i)) * Dx;
}

// Getters for cell volumes

double getCellVol(int i, int j) {
	return cellVol[i][j];
}

// This getter function is NOT compatible with flat control volumes
double getCellVol(int i) {
	return cellVol[i][1];
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Cell Side Normal Vectors
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Stores the normal vectors of on cell top sides
double cellTopSideNormX[Nx + 1][Np + 1], cellTopSideNormP[Nx + 1][Np + 1];
// We chose to store the top side normal vectors to be consistent with notations in the aritcle

// Calculate the normal vectors on the cell sides
void calcCellSideNorm() {
	for (int i = 1; i <= Nx; ++i)
		for (int j = 0; j <= Np; ++j) {
			double pDiff = getCellTopRightP(i, j) - getCellTopLeftP(i, j);
			double topSideVecSizeInv = 1 / sqrt(pDiff * pDiff + Dx * Dx);
			cellTopSideNormX[i][j] = -pDiff * topSideVecSizeInv;
			cellTopSideNormP[i][j] = Dx * topSideVecSizeInv;
			printf("(%d, %d): %f\n", i, j, cellTopSideNormP[i][j]);
		}
}

// Getters for normal vectors on cell TOP sides. Getters for the bottom sides are not defined
// on purpose, to be consistent with the notations in the article

double getCellTopSideNormX(int i, int j) {
	return cellTopSideNormX[i][j];
}

double getCellTopSideNormP(int i, int j) {
	return cellTopSideNormP[i][j];
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Dp Values at Cell Centers
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

double cellCenterDpVec[numCellsX];  // Stores the Dp values at cell centers x

// Calculate Dp values at cell centers
void calcCellCenterDp() {
	for (int i = 0; i < numCellsX; ++i) {
		double x = getCellCenterX(i);
		cellCenterDpVec[i] = ((*pB_fcnPtr)(x) - pA) * NpInv;
	}
}

// Getters for Dp values at cell centers

double getCellCenterDp(int i, int j) {
	return cellCenterDpVec[i];
}

double getCellCenterDp(int i) {
	return cellCenterDpVec[i];
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Set All Mesh Parameters and Build Mesh
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void setMesh() {
	clock_t start = clock();

	calcGridPts();
	calcBaryCenters();
	calcCellSideLen();
	calcCellVol();
	calcCellSideNorm();
	calcCellCenterDp();

	printf("- Mesh setting complete. Time used = %1.2fs.\n",
			((double) (clock() - start)) / CLOCKS_PER_SEC);
}

#endif /* MESH_H_ */
