/*
 * MathFcns.h
 *
 *  Created on: Jun 10, 2017
 *      Author: chuckjia
 *
 *  This file contains mathematical functions that are independent from
 *  model settings.
 */

#ifndef MATHFCNS_H_
#define MATHFCNS_H_
#include "TestFcns.h"

/*
 * The method realizes the mathematical function min(x, y)
 */
double min(double x, double y) {
	if (x < y)
		return x;
	return y;
}

/*
 * The method realizes the mathematical function max(x, y)
 */
double max(double x, double y) {
	if (x > y)
		return x;
	return y;
}

/*
 * The method realizes the mathematical function min(x, y, z)
 */
double min3(double x, double y, double z) {
	if (x < y) {
		if (x < z)
			return x;
		return z;
	}
	// when y <= x
	if (y < z)
		return y;
	return z;
}

/*
 * The method realizes the mathematical function max(x, y, z)
 */
double max3(double x, double y, double z) {
	if (x > y) {
		if (x > z)
			return x;
		return z;
	}
	// when y >= x
	if (y > z)
		return y;
	return z;
}

/*
 * The method realizes the mathematical function minmod(x, y, z)
 */
double minmod3(double x, double y, double z) {
	if (x > 0 && y > 0 && z > 0)
		return min3(x, y, z);
	else if (x < 0 && y < 0 && z < 0)
		return max3(x, y, z);
	else
		return 0;
}

/*
 * The method realizes the sign function sgn(x, y, z)
 */
double sign(double x) {
	if (x > 0)
		return 1;
	if (x < 0)
		return -1;
	return 0;
}

#endif /* MATHFCNS_H_ */
