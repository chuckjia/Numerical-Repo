/*
 * MeshMatrix.h
 *
 *  Created on: Mar 12, 2018
 *      Author: chuckjia
 */

#ifndef INCLUDE_MESHGRID_H_
#define INCLUDE_MESHGRID_H_

#include "Matrix2D.h"

// This class stores the x or y coordinates of the grid point at the cell bottom left corners.
class MeshGrid : public Matrix2D {
public:

	MeshGrid() : Matrix2D() { }

	MeshGrid(int _nrow, int _ncol) : Matrix2D(_nrow, _ncol) { }

	/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
	 * Informational
	 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

	inline double bottLeft(int i, int j) const { return data_slices_[i][j]; }

	inline double bottRight(int i, int j) const { return data_slices_[i + 1][j]; }

	inline double topLeft(int i, int j) const { return data_slices_[i][j + 1]; }

	inline double topRight(int i, int j) const { return data_slices_[i + 1][j + 1]; }
};

#endif /* INCLUDE_MESHGRID_H_ */
