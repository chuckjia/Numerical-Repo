/*
 * Models.h
 *
 *  Created on: Jul 27, 2017
 *      Author: chuckjia
 */

#ifndef MODELS_H_
#define MODELS_H_
#include "Mesh.h"

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Model Functions
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

const double c1_uFcn = M_PI / pB, c2_uFcn = 2 * M_PI / xf;

double u_fcn(double x, double p) {
	return 7.5 + cos(c1_uFcn * p) * cos(c2_uFcn * x);
}

double u_xDer_fcn(double x, double p) {
	return -cos(c1_uFcn * p) * c2_uFcn * sin(c2_uFcn * x);
}

const double c1_wFcn = 2 * M_PI / pB, c2_wFcn = 2 * M_PI / xf;

double w_fcn(double x, double p) {
	return sin(c1_wFcn * p) * cos(c2_wFcn * x);
}

double w_pDer_fcn(double x, double p) {
	return c1_wFcn * cos(c1_wFcn * p) * cos(c2_wFcn * x);
}

const double c1_Test1 = 0.2 * 2 * M_PI;
const double c2_Test1 = 0.2 * 2 * M_PI;
const double c3_Test1 = 2 * M_PI;

double exact_T_Test1(double x, double p, double t) {
	return sin(c1_Test1 * x) * sin(c2_Test1 * p) * cos(c3_Test1 * t);
}

double boundaryVal = 0;

double source_Test1(double T, double x, double p, double t) {
	double xInput = c1_Test1 * x, pInput = c2_Test1 * p, tInput = c3_Test1 * t;
	double xPart = sin(xInput), pPart = sin(pInput), tPart = cos(tInput);
	return - xPart * pPart * c3_Test1 * sin(tInput) +
			T * (u_xDer_fcn(x, p) + w_pDer_fcn(x, p)) +
			c1_Test1 * cos(c1_Test1 * x) * pPart * tPart * u_fcn(x, p) +
			xPart * c2_Test1 * cos(c2_Test1 * p) * tPart * w_fcn(x, p);
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Initial Conditions
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

double soln[numCellsX][numCellsP];

void enforceInitCond() {
	for (int i = 0; i < numCellsX; i++)
		for (int j = 0; j < numCellsP; j++) {
			double x = getCellCenterX(i, j), p = getCellCenterP(i, j);
			soln[i][j] = exact_T_Test1(x, p, 0);
		}
}

#endif /* MODELS_H_ */
