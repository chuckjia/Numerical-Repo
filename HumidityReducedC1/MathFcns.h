/*
 * MathFcns.h
 *
 *  Created on: Jun 10, 2017
 *      Author: chuckjia
 */

#ifndef MATHFCNS_H_
#define MATHFCNS_H_
#include "TestFcns.h"

double min3(double x, double y, double z) {
	if (x < y) {
		if (x < z)
			return x;
		return z;
	} else { // when y <= x
		if (y < z)
			return y;
		return z;
	}
}

double max3(double x, double y, double z) {
	if (x > y) {
		if (x > z)
			return x;
		return z;
	} else { // when y >= x
		if (y > z)
			return y;
		return z;
	}
}

double minmod3(double x, double y, double z) {
	if (x > 0 && y > 0 && z > 0)
		return min3(x, y, z);
	else if (x < 0 && y < 0 && z < 0)
		return max3(x, y, z);
	else
		return 0;
}

double sign(double x) {
	if (x > 0)
		return 1;
	if (x < 0)
		return -1;
	return 0;
}

#endif /* MATHFCNS_H_ */
