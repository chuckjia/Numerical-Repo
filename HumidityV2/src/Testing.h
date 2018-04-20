/*
 * Testing.h
 *
 *  Created on: Mar 13, 2018
 *      Author: chuckjia
 */

#ifndef SRC_TESTING_H_
#define SRC_TESTING_H_

#include "../include/Vector1D.h"

void test() {
	double arr[] = {1.1, 2,3,4,5,6,7,8,9,10};
	Vector1D a(10), b(20);
	for (int i = 0; i < 10; ++i)
		a[i] = arr[i];
	b.print_to_console();
	b = a;
	b.print_to_console();
}

#endif /* SRC_TESTING_H_ */
