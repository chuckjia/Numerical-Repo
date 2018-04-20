/*
 * Vector1D.h
 *
 *  Created on: Apr 2, 2018
 *      Author: chuckjia
 */

#ifndef INCLUDE_VECTOR1D_H_
#define INCLUDE_VECTOR1D_H_

#include <string>
#include "OneDimArray.h"


class Vector1D : public OneDimArray<double> {
public:
	Vector1D() : OneDimArray<double>() { }

	Vector1D(int _len) : OneDimArray<double>(_len) {
		memset(dataPtr(), 0, sizeof(double) * _len);
	}

	Vector1D(int _len, double* arr) : OneDimArray<double>(_len, arr) { }

	void print_to_console() {
		int last = len() - 1;
		printf("[ ");
		for (int i = 0; i < last; ++i)
			printf("%1.4e, ", data_area_[i]);
		printf("%1.4e ]\n", data_area_[last]);
	}
};

#endif /* INCLUDE_VECTOR1D_H_ */
