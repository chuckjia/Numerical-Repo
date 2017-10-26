/*
 * Mesh.h
 *
 *  Created on: Oct 11, 2017
 *      Author: chuckjia
 */

#ifndef MESH_H_
#define MESH_H_
#include "Models.h"

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Auto-generated Constants
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

const int Nx = numDivisions;  // Number of divisions in x direction in space
const int Np = numDivisions;  // Number of divisions in p direction in space

double NxInv = 1. / Nx, NpInv = 1. / Np;

// Following constants are set up for code readability
const int numCellsX = Nx + 2;  // Number of cells in the x direction
const int numCellsP = Np + 2;  // Number of cells in the p direction
int lastRealIndexX = Nx;  // x direction index of cells on the right side of domain
int lastRealIndexP = Np;  // p direction index of cells on the top side of domain
int lastGhostIndexX = Nx + 1;  // x direction index of ghost cells on right side of domain
int lastGhostIndexP = Np + 1;  // p direction index of ghost cells on top side of domain
const int numGridPtsX = numCellsX + 1;  // Number of grid points in x direction (including ghosts)
const int numGridPtsP = numCellsP + 1;  // Number of grid points in p direction (including ghosts)
const int NxPlusOne = Nx + 1;
const int NpPlusOne = Np + 1;

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Grid Points
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

double Dx, DxInv;  // x step size, defined in calcGridPts()
double cellLeftDpVec[numGridPtsX];
double meshgridX[numGridPtsX];
double meshgridP[numGridPtsX][numGridPtsP];

void calcGridPts() {
	Dx = (xf - x0) / Nx;
	DxInv = Nx / (xf - x0);
	for (int i = 0; i < numGridPtsX; ++i) {
		double x = x0 + (i - 1) * Dx;
		meshgridX[i] = x;
		double Dp = ((*pB_fcnPtr)(x) - pA) / Np;
		cellLeftDpVec[i] = Dp;
		for (int j = 0; j < numGridPtsP; ++j)
			meshgridP[i][j] = pA + (j - 1) * Dp;
	}
}

double getCellBottSideLen(int i, int j) {
	return Dx;
}

double getCellBottSideLen() {
	return Dx;
}

double getCellTopSideLen(int i, int j) {
	return Dx;
}

double getCellTopSideLen() {
	return Dx;
}

double getCellLeftDp(int i, int j) {
	return cellLeftDpVec[i];
}

double getCellLeftDp(int i) {
	return cellLeftDpVec[i];
}

double getCellLeftSideLen(int i, int j) {
	return cellLeftDpVec[i];
}

double getCellLeftSideLen(int i) {
	return cellLeftDpVec[i];
}

double getCellRightDp(int i, int j) {
	return cellLeftDpVec[i + 1];
}

double getCellRightDp(int i) {
	return cellLeftDpVec[i + 1];
}

double getCellRightSideLen(int i, int j) {
	return cellLeftDpVec[i + 1];
}

double getCellRightSideLen(int i) {
	return cellLeftDpVec[i + 1];
}

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

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Cell Center
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

double cellCenterX[numCellsX], cellCenterP[numCellsX][numCellsP];

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

double getCellCenterX(int i, int j) {
	return cellCenterX[i];
}

double getCellCenterX(int i) {
	return cellCenterX[i];
}

double getCellCenterP(int i, int j) {
	return cellCenterP[i][j];
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Cell Volumes
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

double cellVol[numCellsX][numCellsP];  // Made 2D to accommodate flat control volumes

void calcCellVol() {
	for (int i = 0; i < numCellsX; ++i)
		for (int j = 0; j < numCellsP; ++j)
			cellVol[i][j] = 0.5 * (getCellLeftDp(i) + getCellRightDp(i)) * Dx;
}

double getCellVol(int i, int j) {
	return cellVol[i][j];
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Cell Side Normal Vectors
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

double cellTopSideNormX[Nx + 1][Np + 1], cellTopSideNormP[Nx + 1][Np + 1];

void calcCellSideNorm() {
	for (int i = 1; i <= lastRealIndexX; ++i)
		for (int j = 0; j <= lastRealIndexP; ++j) {
			double pDiff = getCellTopRightP(i, j) - getCellTopLeftP(i, j);
			double normSize = sqrt(pDiff * pDiff + Dx * Dx);
			cellTopSideNormX[i][j] = -pDiff / normSize;
			cellTopSideNormP[i][j] = Dx / normSize;
		}
}

double getCellTopSideNormX(int i, int j) {
	return cellTopSideNormX[i][j];
}

double getCellTopSideNormP(int i, int j) {
	return cellTopSideNormP[i][j];
}

double getCellBottSideNormX(int i, int j) {
	return cellTopSideNormX[i][j - 1];
}

double getCellBottSideNormP(int i, int j) {
	return cellTopSideNormP[i][j - 1];
}

double cellCenterDpVec[numCellsX];

void calcCellCenterDp() {
	for (int i = 0; i < numCellsX; ++i) {
		double x = getCellCenterX(i);
		cellCenterDpVec[i] = ((*pB_fcnPtr)(x) - pA) * NpInv;
	}
}

double getCellCenterDp(int i, int j) {
	return cellCenterDpVec[i];
}

double getCellCenterDp(int i) {
	return cellCenterDpVec[i];
}

// Wrapper function to build all mesh values
void setMesh() {
	calcGridPts();
	calcBaryCenters();
	calcCellVol();
	calcCellSideNorm();
	calcCellCenterDp();
}

#endif /* MESH_H_ */
