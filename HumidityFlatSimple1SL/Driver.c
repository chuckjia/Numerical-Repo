/*
 * Driver.c
 *
 *  Created on: Jul 6, 2017
 *      Author: chuckjia
 */

#include "TestFcns.h"

int main() {
	forwardEuler();
	printMsg();
	showL2Errors();
	writeResToFile();
	//printFlux(GG_Top);
}