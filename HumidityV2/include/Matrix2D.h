/*
 * CellCenterValueArray.h
 *
 *  Created on: Mar 5, 2018
 *      Author: chuckjia
 */

#ifndef INCLUDE_MATRIX2D_H_
#define INCLUDE_MATRIX2D_H_

#include <string>
#include "TwoDimArray.h"
#include <math.h>

using namespace std;

// This class represents a 2D double array, which initializes as zeros.
class Matrix2D : public TwoDimArray<double> {
public:
	Matrix2D() : TwoDimArray<double>() { }

	Matrix2D(int _nrow, int _ncol) : TwoDimArray<double> (_nrow, _ncol) {
		memset(dataPtr(), 0, sizeof(double) * _nrow * _ncol);
	}

	// Print matrix to a CSV file
	void printToFile(string filename) {
		int _nrow = nrow(), _ncol = ncol();
		char fullname[50];
		sprintf(fullname, "Output/%s.csv", filename.c_str());

		FILE *f = fopen(fullname, "wb");
		for (int i = 0; i < _nrow; ++i) {
			for (int j = 0; j < _ncol; ++j)
				fprintf(f, "%1.20e,", data_slices_[i][j]);
			fprintf(f, "\n");
		}
		fclose(f);
	}
};

#endif /* INCLUDE_MATRIX2D_H_ */
