/*
 * InitialCondtions.h
 *
 *  Created on: Jun 10, 2017
 *      Author: chuckjia
 */

#ifndef MODELCONDITIONS_H_
#define MODELCONDITIONS_H_
#include "Mesh.h"

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Manufactured solution 2
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

const double suppLeft = x0 + 4. * (xL - x0) / 9.,
		suppRight = x0 + 5. * (xL - x0) / 9.,
		pBValAt0 = 1000,
		pSec = (pBValAt0 - pA) / 9.,
		suppBottom = pA + 4. * pSec,
		suppTop = suppBottom + pSec;

double infSmoothFcn(double x, double x1, double x2) {
    double xTemp = (2. * x - x1 - x2) / (x2 - x1);
    return exp(-1. / (1. - xTemp * xTemp));
}

double infSmoothFcnDer(x, x1, x2) {
    double xTemp = (2. * x - x1 - x2) / (x2 - x1);
    return -4. * (2 * x - x1 - x2) / pow(x2 - x1, 2.)
        * 1. / pow(1. - xTemp * xTemp, 2.) * exp(-1. / (1. - xTemp * xTemp));
}

double TManufactured2(double x, double p, double t) {
	if (x >= suppLeft && x <= suppRight && p >= suppBottom && p <= suppTop) {
		double xTemp = (2. * x - suppLeft - suppRight) / (suppRight - suppLeft);
		double pTemp = (2. * p - suppBottom - suppTop) / (suppTop - suppBottom);
		return exp(t) * exp(-1. / (1. - xTemp * xTemp))
				* exp(-1. / (1. - pTemp * pTemp));
	}
	return 0;
}

double qManufactured2(double x, double p, double t) {
	return 0;
}

void sourceManufactured2(double ans[2], double T, double q, double x, double p, double t) {
	if (x >= suppLeft && x <= suppRight && p >= suppBottom && p <= suppTop) {
		double xTemp = (2. * x - suppLeft - suppRight) / (suppRight - suppLeft);
		double pTemp = (2. * p - suppBottom - suppTop) / (suppTop - suppBottom);
		double part_dTdt = exp(t) * exp(-1. / (1 - xTemp * xTemp)) *
				exp(-1. / (1 - pTemp * pTemp));
		double part_dudx = - cos(M_PI * p / p0_CONST) * 2. * M_PI / xL *
				sin(2. * M_PI * x / xL);
		double part_dTdx = exp(t) * infSmoothFcnDer(x, suppLeft, suppRight) *
				infSmoothFcn(p, suppBottom, suppTop);
		double part_domegadp = M_PI / p0_CONST * cos(M_PI * p / p0_CONST) *
				cos(2. * M_PI * x / xL);
		double part_dTdp = exp(t) * infSmoothFcn(x, suppLeft, suppRight) *
				infSmoothFcnDer(p, suppBottom, suppTop);
		ans[0] = part_dTdt + part_dudx * T + uFcn(x, p) * part_dTdx
				+ part_domegadp * T + omegaFcn(x, p) * part_dTdp;
		ans[1] = 0;
	} else {
		ans[0] = 0;
		ans[1] = 0;
	}
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * The functions from the original model
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

double esFcn(double T) {
	return 6.112 * exp(17.67 * (T - 273.15) / (T - 29.65));
}

void modelSource(double ans[2], double T, double q, double x, double p) {
	double omegaVal = omegaFcn(x, p);
	double qsVal = 0.622 * esFcn(T) / p;
	double deltaVal = 0.25 * (1 - sign(omegaVal)) * (1 + sign(q - qsVal));
	double LVal = 2.5008 * 1e6 - 2.3 * 1000 * (T - 275);
	double FVal = qsVal * T * (LVal * R_CONST - Cp_CONST * Rv_CONST * T) /
			(Cp_CONST * Rv_CONST * T * T + qsVal * LVal * LVal);
	ans[0] = omegaVal / (p * Cp_CONST) * (R_CONST * T - deltaVal * LVal * FVal);
	ans[1] = deltaVal * omegaVal * FVal / p;
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Initial conditions: wrapper
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void setInitCond() {
	for (int i = 0; i < numCellsXDir; i++)
		for (int j = 0; j < numCellsPDir; j++) {
			double xVal = getCenterX(i, j), pVal = getCenterP(i, j);
			soln[i][j][0] = TManufactured2(xVal, pVal, 0);
			soln[i][j][1] = qManufactured2(xVal, pVal, 0);
		}
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Source function: wrapper
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

void sourceFcn(double ans[2], double T, double q, double x, double p, double t) {
	sourceManufactured2(ans, T, q, x, p, t);
}


#endif /* MODELCONDITIONS_H_ */
