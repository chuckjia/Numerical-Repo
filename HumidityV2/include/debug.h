/*
 * debug.h
 *
 *  Created on: Mar 14, 2018
 *      Author: chuckjia
 */

#ifndef INCLUDE_DEBUG_H_
#define INCLUDE_DEBUG_H_

#include <stdio.h>
#include <string>

using namespace std;

void tm() {
	printf("\n===== ===== ===== ===== ===== ===== \n");
	printf(">> The program passed here.\n");
	printf("===== ===== ===== ===== ===== ===== \n\n");
}

void tm(string msg) {
	printf("\n===== ===== ===== ===== ===== ===== \n");
	printf(">> The program passed here.\n");
	printf(">> Message: %s\n", msg.c_str());
	printf("===== ===== ===== ===== ===== ===== \n\n");
}

void show_int(int i) {
	printf("\nThe integer value = %d\n", i);
}

void show_int(int i, int j) {
	printf("\nThe integer values = (%d, %d)\n", i, j);
}

void show_int(string s, int i) {
	printf("\n%s = %d\n", s.c_str(), i);
}

void show_int(string s1, string s2, int i, int j) {
	printf("\n(%s, %s) = (%d, %d)\n", s1.c_str(), s2.c_str(), i, j);
}


#endif /* INCLUDE_DEBUG_H_ */
