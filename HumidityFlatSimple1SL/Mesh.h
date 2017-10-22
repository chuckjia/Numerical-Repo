/*
 * Mesh.h
 *
 *  Created on: Jul 27, 2017
 *      Author: chuckjia
 */

#ifndef MESH_H_
#define MESH_H_
#include "Constants.h"

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Getter Functions For the Mesh
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

// Return the x coordinate of cell left side
double getCellLeftX(int i, int j) {
	return x0 + (i - 1) * Dx;
}

// Return the x coordinate of cell right side
double getCellRightX(int i, int j) {
	return x0 + i * Dx;
}

// Return the p coordinate of cell bottom side
double getCellBottP(int i, int j) {
	return pA + (j - 1) * Dp;
}

// Return the p coordinate of cell top side
double getCellTopP(int i, int j) {
	return pA + j * Dp;
}

// Return the x coordinate of cell center
double getCellCenterX(int i, int j) {
	return x0 + (i - 0.5) * Dx;
}

// Return the p coordinate of cell center
double getCellCenterP(int i, int j) {
	return pA + (j - 0.5) * Dp;
}

#endif /* MESH_H_ */
