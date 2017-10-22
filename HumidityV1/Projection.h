/*
 * Projection.h
 *
 *  Created on: Oct 19, 2017
 *      Author: chuckjia
 */

#ifndef PROJECTION_H_
#define PROJECTION_H_
#include "Godunov.h"

double a_proj_cache[Nx + 1], b_proj_cache[Nx + 1];

void fillCache_proj() {
	for (int i = 0; i <= Nx; ++i) {
		double x = getCellCenterX(i);
		a_proj_cache[i] = (*pBxDer_fcnPtr)(x);
		b_proj_cache[i] = (*pB_fcnPtr)(x);
	}
}



#endif /* PROJECTION_H_ */
