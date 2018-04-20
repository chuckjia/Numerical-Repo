/*
 * Conditions.h
 *
 *  Created on: Mar 20, 2018
 *      Author: chuckjia
 */

#ifndef SRC_CONDITIONS_H_
#define SRC_CONDITIONS_H_
#include "../include/NumSoln.h"

void set_init_cond(NumSoln &sol, Matrix2D &cell_centers_x, Matrix2D &cell_centers_p,
		double (*fptr)(double, double, double)) {
	int nrow = sol.nrow(), ncol = sol.ncol();
	assert(nrow == cell_centers_x.nrow() && nrow == cell_centers_p.nrow() &&
			ncol == cell_centers_x.ncol() && ncol == cell_centers_p.ncol());

	for (int i = 0; i < nrow; ++i)
		for (int j = 0; j < ncol; ++j) {
			double x = cell_centers_x[i][j], p = cell_centers_p[i][j];
			sol[i][j] = (*fptr)(x, p, 0);
		}
}

void set_init_cond(NumSoln &sol1, NumSoln &sol2, NumSoln &sol3,
		Matrix2D &cell_centers_x, Matrix2D &cell_centers_p, double (*fptr)(double, double, double)) {
	set_init_cond(sol1, cell_centers_x, cell_centers_p, fptr);
	set_init_cond(sol2, cell_centers_x, cell_centers_p, fptr);
	set_init_cond(sol3, cell_centers_x, cell_centers_p, fptr);
}


#endif /* SRC_CONDITIONS_H_ */
