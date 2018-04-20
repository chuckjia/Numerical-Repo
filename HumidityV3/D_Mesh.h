/*
 * Mesh.h
 *
 *  Created on: Oct 11, 2017
 *      Author: chuckjia
 */

#ifndef D_MESH_H_
#define D_MESH_H_
#include "C_Models.h"

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Mesh Constants: Auto-generated
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

const int Nx = numDivisions;  // Number of space divisions in the x direction
const int Np = numDivisions;  // Number of space divisions in the p direction

// Following constants are set up for code readability
const int numCellsX = Nx + 2;  // Number of cells in the x direction
const int numCellsP = Np + 2;  // Number of cells in the p direction
int lastRealIndexX = Nx;  // Largest x-index of all non-ghost cells
int lastRealIndexP = Np;  // Largest p-index of all non-ghost cells
int lastGhostIndexX = Nx + 1;  // x-index of ghost cells on the right side of domain
int lastGhostIndexP = Np + 1;  // p-index of ghost cells on the top side of domain
const int numGridPtsX = numCellsX + 1;  // Num of grid points in x direction (including ghosts)
const int numGridPtsP = numCellsP + 1;  // Num of grid points in p direction (including ghosts)

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Numerical Solutions: 2D Arrays
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

double T_sl[numCellsX][numCellsP], q_sl[numCellsX][numCellsP];
double u_sl[numCellsX][numCellsP], w_sl[numCellsX][numCellsP];
double phix_sl[numCellsX][numCellsP];

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Grid Points
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

double Dx;  // x step size, defined in calcGridPts()
double cellLeftDpVec[numGridPtsX];  // Stores the Dp values on cell left sides

// The following two arrays store the grid point coordinates. meshgridX stores x coord on cell
// left sides, and meshgridP stores p coord on cell bottom-left sides
double meshGridX[numGridPtsX];  // Stores x_{i-1/2,j}, independent of j
double meshGridP[numGridPtsX][numGridPtsP];  // Stores p_{i-1/2,j-1/2}

// Calculate and store grid points coordinates in cache
void calcGridPts() {
	Dx = (xf - x0) / Nx;
	for (int i = 0; i < numGridPtsX; ++i) {
		double x = x0 + (i - 1) * Dx;
		meshGridX[i] = x;
		double Dp = ((*pB_fcnPtr)(x) - pA) / Np;
		cellLeftDpVec[i] = Dp;
		for (int j = 0; j < numGridPtsP; ++j)
			meshGridP[i][j] = pA + (j - 1) * Dp;
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
	return meshGridX[i];
}

double getCellLeftX(int i) {
	return meshGridX[i];
}

double getCellRightX(int i, int j) {
	return meshGridX[i + 1];
}

double getCellRightX(int i) {
	return meshGridX[i + 1];
}

// Getters for p coordinates on grid points

double getCellBottLeftP(int i, int j) {
	return meshGridP[i][j];
}

double getCellBottRightP(int i, int j) {
	return meshGridP[i + 1][j];
}

double getCellTopLeftP(int i, int j) {
	return meshGridP[i][j + 1];
}

double getCellTopRightP(int i, int j) {
	return meshGridP[i + 1][j + 1];
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Cell Centers
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Store the coordinates for cell barycenters
double cellCenterX[numCellsX], cellCenterP[numCellsX][numCellsP];

// Calculate the coordinates of the intersection of 2 lines. For each line, the coordinates of two points are given.
// Assume (x1, p1) and (x3, p3) are on line 1 with slope k1; (x2, p2) and (x4, p4) are on line 2 with slope k2. See details in notes
// Calculates the coordinates of the intersection point (x, y) = (ans[0], ans[1])
void calcLineIntersection(double ans[2], double x1, double x2, double x3, double x4, double p1, double p2, double p3, double p4) {
	double k1 = (p3 - p1) / (x3 - x1), k2 = (p4 - p2) / (x4 - x2);

	double x = (k1 * x1 - p1 - k2 * x2 + p2) / (k1 - k2);
	ans[0] = x;
	ans[1] = k1 * (x - x1) + p1;
}

// Assume that (x1, p1), (x2, p2), (x3, p3), and (x4, p4) are arranged in counter-clockwise order
void calcTrapezoidCenter(double center[2], double x1, double x2, double x3, double x4, double p1, double p2, double p3, double p4) {
	double x124 = (x1 + x2 + x4) / 3, x134 = (x1 + x3 + x4) / 3, x123 = (x1 + x2 + x3) / 3, x234 = (x2 + x3 + x4) / 3;
	double p124 = (p1 + p2 + p4) / 3, p134 = (p1 + p3 + p4) / 3, p123 = (p1 + p2 + p3) / 3, p234 = (p2 + p3 + p4) / 3;

	calcLineIntersection(center, x124, x123, x234, x134, p124, p123, p234, p134);
}

// Calculate cell barycenters. Require cells to be trapezoids
void calcBaryCenters() {
	for (int i = 0; i < numCellsX; ++i)
		for (int j = 0; j < numCellsP; ++j) {
			double center[2];
			double x1 = getCellLeftX(i, j), x2 = getCellRightX(i, j), x3 = x2, x4 = x1;
			double p1 = getCellBottLeftP(i, j), p2 = getCellBottRightP(i, j), p3 = getCellTopRightP(i, j), p4 = getCellTopLeftP(i, j);
			calcTrapezoidCenter(center, x1, x2, x3, x4, p1, p2, p3, p4);
			cellCenterX[i] = center[0];
			cellCenterP[i][j] = center[1];
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
	for (int i = 0; i < numCellsX; ++i)
		for (int j = 0; j < numCellsP; ++j)
			cellBottSideLenVec[i][j] = sqrt(Dx * Dx + pow(getCellBottRightP(i, j) - getCellBottLeftP(i, j), 2));
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
			double topSideVecSize = sqrt(pDiff * pDiff + Dx * Dx);
			cellTopSideNormX[i][j] = -pDiff / topSideVecSize;
			cellTopSideNormP[i][j] = Dx / topSideVecSize;
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
		cellCenterDpVec[i] = ((*pB_fcnPtr)(x) - pA) /Np;
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

#endif /* D_MESH_H_ */
