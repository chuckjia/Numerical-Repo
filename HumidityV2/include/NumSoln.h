/*
 * NumSolnMatrix.h
 *
 *  Created on: Mar 13, 2018
 *      Author: chuckjia
 */

#ifndef INCLUDE_NUMSOLN_H_
#define INCLUDE_NUMSOLN_H_

#include "Matrix2D.h"
#include "vector"

class NumSoln : public Matrix2D {
public:
	NumSoln() : Matrix2D() { }

	NumSoln(int _nrow, int _ncol) : Matrix2D(_nrow, _ncol) { }

	void enforce_dirichelt_left(double bd_val) {
		int n = ncol_;
		double twice_bd_val = 2 * bd_val;
		for (int j = 0; j < n; ++j)
			data_slices_[0][j] = twice_bd_val - data_slices_[1][j];
	}

	void enforce_dirichelt_left(vector<double> &bd_val) {
		int n = ncol_;
		for (int j = 0; j < n; ++j)
			data_slices_[0][j] = 2 * bd_val[j] - data_slices_[1][j];
	}

	void enforce_nuemann_right() {
		int n = ncol_;
		for (int j = 0; j < n; ++j)
			data_slices_[0][j] = data_slices_[1][j];
	}
};

#endif /* INCLUDE_NUMSOLN_H_ */
