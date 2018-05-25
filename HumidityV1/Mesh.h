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

const int Nx = numDivision;  // Number of space divisions in the x direction
const int Np = numDivision;  // Number of space divisions in the p direction

// Following constants are set up for code readability
double NxInv = 1. / Nx, NpInv = 1. / Np;
const int numCellX = Nx + 2;  // Number of cells in the x direction
const int numCellP = Np + 2;  // Number of cells in the p direction
int lastRealIndexX = Nx;  // Largest x-index of all non-ghost cells
int lastRealIndexP = Np;  // Largest p-index of all non-ghost cells
int lastGhostIndexX = Nx + 1;  // x-index of ghost cells on the right side of domain
int lastGhostIndexP = Np + 1;  // p-index of ghost cells on the top side of domain
const int numGridPtX = numCellX + 1;  // Num of grid points in x direction (including ghosts)
const int numGridPtP = numCellP + 1;  // Num of grid points in p direction (including ghosts)
const int NxPlusOne = Nx + 1;
const int NpPlusOne = Np + 1;

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Grid Points
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

double Dx, DxInv;  // x step size, defined in calcGridPts()
double cellLeftDp_[numGridPtX];  // Stores the Dp values on cell left sides

// The following two arrays store the grid point coordinates. meshgridX stores x coord on cell
// left sides, and meshgridP stores p coord on cell bottom-left sides
double meshGridX_[numGridPtX];  // Stores x_{i-1/2,j}
double meshGridP_[numGridPtX][numGridPtP];  // Stores p_{i-1/2,j-1/2}

// Calculate and store grid points coordinates in cache
void calcGridPts() {
	Dx = (xf - x0) / Nx;
	DxInv = Nx / (xf - x0);
	for (int i = 0; i < numGridPtX; ++i) {
		double x = x0 + (i - 1) * Dx;
		meshGridX_[i] = x;
		double Dp = ((*pB_fptr)(x) - pA) * NpInv;
		cellLeftDp_[i] = Dp;
		for (int j = 0; j < numGridPtP; ++j)
			meshGridP_[i][j] = pA + (j - 1) * Dp;
	}
}

// Getters for the Dp values

double getCellLeftDp(int i, int j) {
	return cellLeftDp_[i];
}

double getCellLeftDp(int i) {
	return cellLeftDp_[i];
}

double getCellRightDp(int i, int j) {
	return cellLeftDp_[i + 1];
}

double getCellRightDp(int i) {
	return cellLeftDp_[i + 1];
}

// Getters for x coordinates on grid points

double getCellLeftX(int i, int j) {
	return meshGridX_[i];
}

double getCellLeftX(int i) {
	return meshGridX_[i];
}

double getCellRightX(int i, int j) {
	return meshGridX_[i + 1];
}

double getCellRightX(int i) {
	return meshGridX_[i + 1];
}

// Getters for p coordinates on grid points

double getCellBottLeftP(int i, int j) {
	return meshGridP_[i][j];
}

double getCellBottRightP(int i, int j) {
	return meshGridP_[i + 1][j];
}

double getCellTopLeftP(int i, int j) {
	return meshGridP_[i][j + 1];
}

double getCellTopRightP(int i, int j) {
	return meshGridP_[i + 1][j + 1];
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Cell Centers
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Store the coordinates for cell barycenters
double cellCenterX_[numCellX], cellCenterP_[numCellX][numCellP];

// Calculate cell barycenters. Require cells to be trapezoids
void calcBaryCenters() {
	double h = Dx / 3.;  // h = (x2_Outer - x1_Outer) / 3
	for (int i = 0; i < numCellX; ++i) {
		double x1_Outer = getCellLeftX(i), x2_Outer = getCellRightX(i);
		double x1_Inner = (2 * x1_Outer + x2_Outer) * ONE_THIRD;
		double a = getCellRightDp(i) * ONE_THIRD,  // a = (p3_Outer - p2_Outer) / 3
				b = getCellLeftDp(i) * ONE_THIRD,  // b = (p4_Outer - p1_Outer) / 3
				inv_ab_sum = 1. / (a + b);
		cellCenterX_[i] = x1_Inner + h * a * inv_ab_sum;
		for (int j = 0; j < numCellP; ++j) {
			double p1_Outer = getCellBottLeftP(i, j), p2_Outer = getCellBottRightP(i, j),
					p3_Outer = getCellTopRightP(i, j), p4_Outer = getCellTopLeftP(i, j);
			double p1_Inner = (p1_Outer + p2_Outer + p4_Outer) * ONE_THIRD;
			double bLeft = (p4_Outer - p3_Outer) * ONE_THIRD;
			cellCenterP_[i][j] = p1_Inner + (b - bLeft) * a * inv_ab_sum;
		}
	}
}

// Getters for coordinates of cell centers

double getCellCenterX(int i, int j) {
	return cellCenterX_[i];
}

double getCellCenterX(int i) {
	return cellCenterX_[i];
}

double getCellCenterP(int i, int j) {
	return cellCenterP_[i][j];
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Cell Side Lengths
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Stores the length of cell bottom sides
double cellBottSideLen_[numCellX][numCellP];

// Calculate all cell sides length: calculate bottom side length of cell (i, j)
void calcCellSideLen() {
	for (int i = 0; i < numCellX; ++i) {
		double xDiffSq = pow(getCellRightX(i) - getCellLeftX(i), 2);
		for (int j = 0; j < numCellP; ++j) {
			cellBottSideLen_[i][j] = sqrt(xDiffSq +
					pow(getCellBottRightP(i, j) - getCellBottLeftP(i, j), 2));
		}
	}
}

// Getters for cell side lengths

double getCellBottSideLen(int i, int j) {
	return cellBottSideLen_[i][j];
}

double getCellTopSideLen(int i, int j) {
	return cellBottSideLen_[i][j + 1];
}

double getCellLeftSideLen(int i, int j) {
	return cellLeftDp_[i];
}

double getCellLeftSideLen(int i) {
	return cellLeftDp_[i];
}

double getCellRightSideLen(int i, int j) {
	return cellLeftDp_[i + 1];
}

double getCellRightSideLen(int i) {
	return cellLeftDp_[i + 1];
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Cell Volumes
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Store the volumes of all cells. Made 2D to accommodate flat control volumes
double cellVol_[numCellX][numCellP];

// Calculate cell volumes. Require cells to be trapezoids
void calcCellVol() {
	for (int i = 0; i < numCellX; ++i)
		for (int j = 0; j < numCellP; ++j)
			cellVol_[i][j] = 0.5 * (getCellLeftDp(i) + getCellRightDp(i)) * Dx;
}

// Getters for cell volumes

double getCellVol(int i, int j) {
	return cellVol_[i][j];
}

// This getter function is NOT compatible with flat control volumes
double getCellVol(int i) {
	return cellVol_[i][1];
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Cell Side Normal Vectors
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Stores the normal vectors of on cell top sides
double cellTopSideNormVecX_[Nx + 1][Np + 1], cellTopSideNormVecP_[Nx + 1][Np + 1];
// We chose to store the top side normal vectors to be consistent with notations in the aritcle

// Calculate the normal vectors on the cell sides
void calcCellSideNormVec() {
	for (int i = 1; i <= Nx; ++i)
		for (int j = 0; j <= Np; ++j) {
			double pDiff = getCellTopRightP(i, j) - getCellTopLeftP(i, j);
			double topSideVecSizeInv = 1 / sqrt(pDiff * pDiff + Dx * Dx);
			cellTopSideNormVecX_[i][j] = -pDiff * topSideVecSizeInv;
			cellTopSideNormVecP_[i][j] = Dx * topSideVecSizeInv;
		}
}

// Getters for normal vectors on cell TOP sides. Getters for the bottom sides are not defined
// on purpose, to be consistent with the notations in the article

double getCellTopSideNormVecX(int i, int j) {
	return cellTopSideNormVecX_[i][j];
}

double getCellTopSideNormVecP(int i, int j) {
	return cellTopSideNormVecP_[i][j];
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Dp Values at Cell Centers
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

double cellCenterDp_[numCellX];  // Stores the Dp values at cell centers x

// Calculate Dp values at cell centers
void calcCellCenterDp() {
	for (int i = 0; i < numCellX; ++i) {
		double x = getCellCenterX(i);
		cellCenterDp_[i] = ((*pB_fptr)(x) - pA) * NpInv;
	}
}

// Getters for Dp values at cell centers

double getCellCenterDp(int i, int j) {
	return cellCenterDp_[i];
}

double getCellCenterDp(int i) {
	return cellCenterDp_[i];
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Numerical Solutions: 2D Arrays
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

double T_[numCellX][numCellP], q_[numCellX][numCellP];
double u_[numCellX][numCellP], w_[numCellX][numCellP];
double phix_[numCellX][numCellP];

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Set All Mesh Parameters and Build Mesh
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void setMesh() {
	clock_t start = clock();

	calcGridPts();
	calcBaryCenters();
	calcCellSideLen();
	calcCellVol();
	calcCellSideNormVec();
	calcCellCenterDp();

	printf("- Mesh setting complete. Time used = %1.2fs.\n",
			((double) (clock() - start)) / CLOCKS_PER_SEC);
}

#endif /* MESH_H_ */
