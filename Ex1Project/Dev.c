#include <stdio.h>
#include "probset.h"

int showStep = 1;
int maxItr = 10000;

/* ===== ===== ===== ===== */

void secantAitkensDelta2() {
	printf("- Secant Method With Aitken's Delta 2\n");
	double x0, x1;
	setInit(&x0, &x1, probNum);
	double x2 = (x0 * func(x1) - x1 * func(x0)) / (func(x1) - func(x0));
	double xHat = x0 - pow(x1 - x0, 2) / (x2 - 2 * x1 + x0);

	double fx0Val = func(x0), fx1Val = func(x1), fx2Val = func(x2);
	int k = 0;
	while (k < maxItr && fabs(func(xHat)) > err) {
		k++;
		double denom = fx2Val - fx1Val;
		if (fabs(denom) < err) {
			printf("Here");
			break;
		}
		double xNew = (x1 * fx2Val - x2 * fx1Val) / denom;
		x0 = x1;
		x1 = x2;
		x2 = xNew;
		xHat = x0 - pow(x1 - x0, 2) / (x2 - 2 * x1 + x0);
		fx0Val = func(x0);
		fx1Val = func(x1);
		fx2Val = func(x2);
		if (showStep)
			printf("Iteration %d, Result: %10.13f\n", k, xHat);
	}
	printf("Number of iterations: %d\nResult: %10.13f\n\n", k, xHat);
}

void bisection() {
	printf("- Bisection Method\n");
	double x0, x1;
	setInit(&x0, &x1, probNum);

	int k = 0;
	double fx1Val = func(x1);
	while(k < maxItr && fabs(x1 - x0) > err) {
		k++;
		double xNew = (x0 + x1) * 0.5;
		double fxNewVal = func(xNew);
		if (fx1Val * fxNewVal < 0) {
			x0 = xNew;
		} else {
			x1 = xNew;
		}
		fx1Val = func(x1);
		if (showStep)
			printf("Iteration %d, Result: %10.13f, %10.13f\n", k, x0, x1);
	}
	printf("Number of iterations: %d\nResult: %10.13f\n\n", k, (x0 + x1) * 0.5);
}

void IQI() {
	printf("- Inverse Quadratic Interpolation Method\n");
	double x0, x1;
	setInit(&x0, &x1, probNum);
	double x2 = (x0 * func(x1) - x1 * func(x0)) / (func(x1) - func(x0));

	int k = 0;
	double fx0Val = func(x0), fx1Val = func(x1), fx2Val = func(x2);
	while(k < maxItr && fabs(fx2Val) > err) {
		k++;
		double xNew;
		double denom = (fx0Val - fx1Val) * (fx0Val - fx2Val);
		if (fabs(denom) < err) break;
		xNew = fx1Val * fx2Val / denom * x0;
		denom = (fx1Val - fx0Val) * (fx1Val - fx2Val);
		if (fabs(denom) < err) break;
		xNew += fx0Val * fx2Val / denom * x1;
		denom = (fx2Val - fx0Val) * (fx2Val - fx1Val);
		if (fabs(denom) < err) break;
		xNew += fx0Val * fx1Val / denom * x2;
		x0 = x1;
		x1 = x2;
		x2 = xNew;
		fx0Val = func(x0);
		fx1Val = func(x1);
		fx2Val = func(x2);
		if (showStep)
			printf("Iteration %d, Result: %10.13f\n", k, x2);
	}
	printf("Number of iterations: %d\nResult: %10.13f\n\n", k, x2);
}

int main() {
	selectProb();
	initFunc();
	printf("===== ===== ===== ===== ===== ===== =====");
	printf("\nProblem %d\n\n", probNum);
	IQI();
}
