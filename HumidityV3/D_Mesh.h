/*
 * D_Mesh.h
 *
 *  Created on: Oct 11, 2017
 *      Author: chuckjia
 *
 *  This file contains global variables and functions that store and build mesh of the physical domain.
 */

#ifndef D_MESH_H_
#define D_MESH_H_
#include "C_Models.h"

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Mesh Constants: Auto-generated
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

const int Nx = numDivision;  // Number of space divisions, x-direction
const int Np = numDivision;  // Number of space divisions, p-direction

// Following constants are set up for code readability
const int numCellX = Nx + 2;  // Number of cells, x-direction
const int numCellP = Np + 2;  // Number of cells, p-direction
int lastRealIndexX = Nx;  // Largest x-index of all non-ghost cells
int lastRealIndexP = Np;  // Largest p-index of all non-ghost cells
int lastGhostIndexX = Nx + 1;  // x-index of ghost cells on the RIGHT side of domain
int lastGhostIndexP = Np + 1;  // p-index of ghost cells on the TOP side of domain
const int numGridPtX = numCellX + 1;  // Number of grid points, x-direction (including ghost cells)
const int numGridPtP = numCellP + 1;  // Number of grid points, p-direction (including ghost cells)


/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Numerical Solutions: 2D Arrays
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

double T_[numCellX][numCellP], q_[numCellX][numCellP];
double u_[numCellX][numCellP], w_[numCellX][numCellP], phix_[numCellX][numCellP];


/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Grid Points
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

double Dx;  // x step size, defined in calcGridPts()
double cellLeftDp_[numGridPtX];  // Stores the Dp values on cell left sides, defined in calcGridPts()

// The following two arrays store the grid point coordinates. meshGridX stores x coord on cell
// left sides, and meshGridP stores p coord on cell bottom-left sides
double meshGridX_[numGridPtX];  // Stores x_{i-1/2,j}, independent of j
double meshGridP_[numGridPtX][numGridPtP];  // Stores p_{i-1/2,j-1/2}

// Calculate and store grid points coordinates in cache
void calcGridPts() {
	Dx = (xf - x0) / Nx;
	for (int i = 0; i < numGridPtX; ++i) {
		double x = x0 + (i - 1) * Dx;
		meshGridX_[i] = x;
		double Dp = ((*pB_fptr)(x) - pA) / Np;
		cellLeftDp_[i] = Dp;
		for (int j = 0; j < numGridPtP; ++j)
			meshGridP_[i][j] = pA + (j - 1) * Dp;
	}
}

/**
 * Getters for the Dp values
 */
double getCellLeftDp(int i, int j) { return cellLeftDp_[i]; }
double getCellLeftDp(int i) { return cellLeftDp_[i]; }

double getCellRightDp(int i, int j) { return cellLeftDp_[i + 1]; }
double getCellRightDp(int i) { return cellLeftDp_[i + 1]; }

/**
 * Getters for x coordinates on grid points
 */
double getCellLeftX(int i, int j) { return meshGridX_[i]; }
double getCellLeftX(int i) { return meshGridX_[i]; }

double getCellRightX(int i, int j) { return meshGridX_[i + 1]; }
double getCellRightX(int i) { return meshGridX_[i + 1]; }

/**
 * Getters for p coordinates on grid points
 */
double getCellBottLeftP(int i, int j) { return meshGridP_[i][j]; }
double getCellBottRightP(int i, int j) { return meshGridP_[i + 1][j]; }
double getCellTopLeftP(int i, int j) { return meshGridP_[i][j + 1]; }
double getCellTopRightP(int i, int j) { return meshGridP_[i + 1][j + 1]; }


/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Cell Centers
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Store the coordinates for cell barycenters
double cellCenterX_[numCellX], cellCenterP_[numCellX][numCellP];

// Finds the intersection of 2 lines. For each line, coordinates of two points are given
// Assume (x1, p1) and (x3, p3) are on line 1 with slope k1; (x2, p2) and (x4, p4) are on line 2 with slope k2. See details in notes
// Calculates the coordinates of the intersection point (x, y) = (ans[0], ans[1])
void calcLineIntersection(double ans[2], double x1, double x2, double x3, double x4, double p1, double p2, double p3, double p4) {
	double k1 = (p3 - p1) / (x3 - x1), k2 = (p4 - p2) / (x4 - x2);
	double x = (k1 * x1 - p1 - k2 * x2 + p2) / (k1 - k2);
	ans[0] = x;
	ans[1] = k1 * (x - x1) + p1;
}

// Calculate the barycenter of each trapezoid cell, using definition in (3.5)
// Assume that (x1, p1), (x2, p2), (x3, p3), and (x4, p4) are arranged in counter-clockwise order
void calcTrapezoidCenter(double center[2], double x1, double x2, double x3, double x4, double p1, double p2, double p3, double p4) {
	double x124 = (x1 + x2 + x4) / 3, x134 = (x1 + x3 + x4) / 3, x123 = (x1 + x2 + x3) / 3, x234 = (x2 + x3 + x4) / 3;
	double p124 = (p1 + p2 + p4) / 3, p134 = (p1 + p3 + p4) / 3, p123 = (p1 + p2 + p3) / 3, p234 = (p2 + p3 + p4) / 3;
	calcLineIntersection(center, x124, x123, x234, x134, p124, p123, p234, p134);
}

// Calculate cell barycenters. Require cells to be trapezoids
void calcBaryCenters() {
	for (int i = 0; i < numCellX; ++i)
		for (int j = 0; j < numCellP; ++j) {
			double center[2];
			double x1 = getCellLeftX(i, j), x2 = getCellRightX(i, j), x3 = x2, x4 = x1;
			double p1 = getCellBottLeftP(i, j), p2 = getCellBottRightP(i, j), p3 = getCellTopRightP(i, j), p4 = getCellTopLeftP(i, j);
			calcTrapezoidCenter(center, x1, x2, x3, x4, p1, p2, p3, p4);
			cellCenterX_[i] = center[0];
			cellCenterP_[i][j] = center[1];
		}
}

/*
 * Getters for coordinates of cell centers
 */
double getCellCenterX(int i, int j) { return cellCenterX_[i]; }
double getCellCenterX(int i) { return cellCenterX_[i]; }

double getCellCenterP(int i, int j) { return cellCenterP_[i][j]; }


/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Cell Side Lengths
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Stores the length of cell bottom sides
double cellBottSideLen_[numCellX][numCellP];

// Calculate all cell sides length: calculate bottom side length of cell (i, j)
void calcCellSideLen() {
	for (int i = 0; i < numCellX; ++i)
		for (int j = 0; j < numCellP; ++j) {
			double diffP = getCellBottRightP(i, j) - getCellBottLeftP(i, j);
			cellBottSideLen_[i][j] = sqrt(Dx * Dx + diffP * diffP);
		}
}

/*
 * Getters for cell side lengths
 */
double getCellBottSideLen(int i, int j) { return cellBottSideLen_[i][j]; }

double getCellTopSideLen(int i, int j) { return cellBottSideLen_[i][j + 1]; }

double getCellLeftSideLen(int i, int j) { return cellLeftDp_[i]; }
double getCellLeftSideLen(int i) { return cellLeftDp_[i]; }

double getCellRightSideLen(int i, int j) { return cellLeftDp_[i + 1]; }
double getCellRightSideLen(int i) { return cellLeftDp_[i + 1]; }


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

/*
 * Getters for cell volumes
 */
double getCellVol(int i, int j) { return cellVol_[i][j]; }

// This getter function is NOT compatible with flat control volumes
double getCellVol(int i) { return cellVol_[i][1]; }


/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Cell Side Normal Vectors
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Stores the normal vectors of on cell top sides
// We chose to store the top side normal vectors to be consistent with notations in the aritcle
double cellTopSideNormVecX_[numCellX][numCellP], cellTopSideNormVecP_[numCellX][numCellP];

// Calculate the normal vectors on the cell sides
void calcCellSideNormVec() {
	for (int i = 1; i <= Nx; ++i)
		for (int j = 0; j <= Np; ++j) {
			double pDiff = getCellTopRightP(i, j) - getCellTopLeftP(i, j),
					topSideNormVecSize = sqrt(pDiff * pDiff + Dx * Dx);
			cellTopSideNormVecX_[i][j] = -pDiff / topSideNormVecSize;
			cellTopSideNormVecP_[i][j] = Dx / topSideNormVecSize;
		}
}

/*
 * Getters for normal vectors on cell TOP sides. Getters for the bottom sides are not defined
 * on purpose, to be consistent with the notations in the article
 */
double getCellTopSideNormVecX(int i, int j) { return cellTopSideNormVecX_[i][j]; }
double getCellTopSideNormVecP(int i, int j) { return cellTopSideNormVecP_[i][j]; }


/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Dp Values at Cell Centers
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

double cellCenterDp_[numCellX];  // Stores the Dp values at cell centers x

// Calculate Dp values at cell centers
void calcCellCenterDp() {
	// Previous approach
	// for (int i = 0; i < numCellX; ++i) { cellCenterDpVec[i] = ((*pB_fptr)(getCellCenterX(i)) - pA) /Np; }
	for (int i = 0; i < numCellX; ++i) {
		double x4 = getCellLeftX(i), x3 = getCellRightX(i),  // Numbering follows the notes
				p4 = getCellTopLeftP(i, Np), p3 = getCellTopRightP(i, Np),
				k = (p4 - p3) / (x4 - x3), xCenter = getCellCenterX(i),
				p = k * (xCenter - x4) + p4;
		cellCenterDp_[i] = (p - pA) / Np;
	}
}

/*
 * Getters for Dp values at cell centers
 */
double getCellCenterDp(int i, int j) { return cellCenterDp_[i]; }
double getCellCenterDp(int i) { return cellCenterDp_[i]; }


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

	printTimeUsed(start, clock(), "ms", "Mesh building completed.");
}

#endif /* D_MESH_H_ */
