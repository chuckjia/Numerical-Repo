/*
 * FV.h
 *
 *  Created on: Oct 8, 2017
 *      Author: chuckjia
 */

#ifndef FV_H_
#define FV_H_
#include "Analysis.h"

void upwind_right() {
	selectModel();
	enforceIC();
	double lambda = Dt / Dx;
	if (printStat)
		printf("t (in sec)	normExact	norm		err		err/norm\n");
	for (int tt = 1; tt <= numTimeSteps; tt++) {
		double t = Dt * tt;
		double norm = 0, err = 0, normExact = 0;
		double uLeft = u[0];
		for (int i = 1; i <= lastRealIndex; i++) {
			double uLeftNext = u[i];
			u[i] += -a * lambda * (u[i] - uLeft);
			uLeft = uLeftNext;
			// Generating statistics
			norm += pow(u[i], 2);
			double x = getCellCenter(i);
			double exactSolnVal = (*exactSolnPtr)(x, t);
			err += pow(u[i] - exactSolnVal, 2);
			normExact += exactSolnVal * exactSolnVal;
		}
		(*enforceBCPtr)();
		if (printStat) {
			norm = sqrt(Dx * norm);
			err = sqrt(Dx * err);
			normExact = sqrt(Dx * normExact);
			printf("t = %2.4fs	%1.6e	%1.6e	%1.6e	%1.6e\n", t, normExact, norm, err, err / norm);
		}
		if (printMovie)
			writeResToFileForMovie(tt);
	}
}

void laxWendroff() {
	selectModel();
	enforceIC();
	double lambda = Dt / Dx;
	if (printStat)
		printf("t (in sec)	normExact	norm		err		err/norm\n");
	for (int tt = 1; tt <= numTimeSteps; tt++) {
		double t = Dt * tt;
		double norm = 0, err = 0, normExact = 0;
		double uLeft = u[0];
		for (int i = 1; i <= lastRealIndex; i++) {
			double uPrev = u[i];
			u[i] += (-a * lambda * (u[i + 1] - uLeft) +
					a * a * lambda * lambda * (u[i + 1] - 2 * u[i] + uLeft)) / 2;
			uLeft = uPrev;
			// Generating statistics
			norm += pow(u[i], 2);
			double x = getCellCenter(i);
			double exactSolnVal = (*exactSolnPtr)(x, t);
			err += pow(u[i] - exactSolnVal, 2);
			normExact += exactSolnVal * exactSolnVal;
		}
		norm = sqrt(Dx * norm);
		err = sqrt(Dx * err);
		normExact = sqrt(Dx * normExact);
		(*enforceBCPtr)();
		if (printStat)
			printf("%t = 2.4fs	%1.6e	%1.6e	%1.6e	%1.6e\n", t, normExact, norm, err, err / norm);
		if (printMovie)
			writeResToFileForMovie(tt);
	}
}

double minmod(double a, double b) {
	if (a * b <= 0)
		return 0;
	if (fabs(a) < fabs(b))
		return a;
	return b;
}

double maxmod(double a, double b) {
	if (a * b <= 0)
		return 0;
	if (fabs(a) > fabs(b))
		return a;
	return b;
}

void laxWendroff_limiter() {
	selectModel();
	enforceIC();
	double lambda = Dt / Dx;
	double v = a * lambda;
	if (printStat)
		printf("t (in sec)	normExact	norm		err		err/norm\n");
	for (int tt = 1; tt <= numTimeSteps; tt++) {
		double t = Dt * tt;
		double norm = 0, err = 0, normExact = 0;
		double uLeft = u[0], u2Left = u[0];
		for (int i = 1; i <= lastRealIndex; i++) {
			double uLeftNext = u[i];
			double slope = maxmod(
					minmod((u[i + 1] - u[i]) / Dx, 2 * (u[i] - uLeft) / Dx),
					minmod(2 * (u[i + 1] - u[i]) / Dx, (u[i] - uLeft) / Dx));
			double slopeLeft = maxmod(
					minmod((u[i] - uLeft) / Dx, 2 * (uLeft - u2Left) / Dx),
					minmod(2 * (u[i] - uLeft) / Dx, (uLeft - u2Left) / Dx)
			);
			u[i] += - a * lambda * (u[i] - uLeft) -
					lambda / 2 * a * (Dx - a * Dt) * (slope - slopeLeft);
			u2Left = uLeft;
			uLeft = uLeftNext;
			// Generating statistics
			norm += pow(u[i], 2);
			double x = getCellCenter(i);
			double exactSolnVal = (*exactSolnPtr)(x, t);
			err += pow(u[i] - exactSolnVal, 2);
			normExact += exactSolnVal * exactSolnVal;
		}
		norm = sqrt(Dx * norm);
		err = sqrt(Dx * err);
		normExact = sqrt(Dx * normExact);
		(*enforceBCPtr)();
		if (printStat)
			printf("t = %2.4fs	%1.6e	%1.6e	%1.6e	%1.6e\n", t, normExact, norm, err, err / norm);
		if (printMovie)
			writeResToFileForMovie(tt);
	}
}

#endif /* FV_H_ */