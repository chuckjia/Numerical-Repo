#include <stdio.h>
#include <time.h>
#include "probset.h"

int showStep = 0;
int maxItr = 10000;

/* ===== ===== ===== ===== ===== ===== ===== =====
 * Newton's method: need derivatives known
 * ===== ===== ===== ===== ===== ===== ===== ===== */

void newton() {
	printf("- Newton's Method\n");
	double x0, x1;
	setInit(&x0, &x1, probNum);

	int k = 0;
	double fx1Val = func(x1);
	while (k < maxItr && fabs(fx1Val) > err) {
		k++;
		double derVal = derFunc(x1);
		if (fabs(derVal) < err)
			break;
		x1 = x1 - fx1Val / derVal;
		fx1Val = func(x1);
		if (showStep)
			printf("Iteration %d, Result: %10.13f\n", k, x1);
	}
	printf("Number of iterations: %d\nResult: %10.13f\n", k, x1);
}

/* ===== ===== ===== ===== ===== ===== ===== =====
 * Secant method
 * ===== ===== ===== ===== ===== ===== ===== ===== */

void secant() {
	printf("- Secant Method\n");
	double x0, x1;
	setInit(&x0, &x1, probNum);

	int k = 0;
	double fx1Val = func(x1);
	while (k < maxItr && fabs(fx1Val) > err) {
		k++;
		double fDer = (fx1Val - func(x0)) / (x1 - x0);
		x0 = x1;
		x1 = x1 - fx1Val / fDer;
		fx1Val = func(x1);
		if (showStep)
			printf("Iteration %d, Result: %10.13f\n", k, x1);
	}
	printf("Number of iterations: %d\nResult: %10.13f\n", k, x1);
}

/* ===== ===== ===== ===== ===== ===== ===== =====
 * Secant method with Regula-Falsi
 * ===== ===== ===== ===== ===== ===== ===== ===== */

void secantRF() {
	printf("- Secant Method With Regula Falsi\n");
	double x0, x1;
	setInit(&x0, &x1, probNum);
	int k = 0;

	double fx0Val = func(x0);
	double fx1Val = func(x1);
	while (k < maxItr && fabs(fx0Val) > err && fabs(fx1Val) > err) {
		k++;
		double denom = fx1Val - fx0Val;
		if (fabs(denom) < err)
			break;
		double xStar = x1 - fx1Val * (x1 - x0) / denom;
		if (func(xStar) * fx0Val < 0) {
			x1 = x0;
			x0 = xStar;
		} else {
			x0 = x1;
			x1 = xStar;
		}
		fx0Val = func(x0);
		fx1Val = func(x1);
		// Print intermediate results
		if (showStep == 1)
			printf("Iteration %d, Result: (%10.13f, %10.13f)\n", k, x0, x1);
	}
	double result = x0;
	if (fabs(fx0Val) > fabs(fx1Val))
		result = x1;
	printf("Number of iterations: %d\nResult: %10.13f\n", k, result);
}

/* ===== ===== ===== ===== ===== ===== ===== =====
 * Ridders' method
 * ===== ===== ===== ===== ===== ===== ===== ===== */

void ridders(){
	printf("- Ridder's Method\n");
	double x0, x2;
	setInit(&x0, &x2, probNum);
	int k = 0;

	double fx0Val = func(x0);
	double fx2Val = func(x2);
	while (k < maxItr && fabs(fx0Val) > err && fabs(fx2Val) > err) {
		k++;
		double x1 = (x0 + x2) / 2;
		double fx1Val = func(x1);
		double x3 = x1 + (x1 - x0) * sgnFunc(fx0Val) * fx1Val / pow(pow(fx1Val,2) - fx0Val * fx2Val, 0.5);
		double fx3Val = func(x3);
		if (fx1Val * fx3Val < 0) {
			x0 = x1;
			x2 = x3;
		} else {
			if (fx0Val * fx3Val < 0)
				x2 = x3;
			else
				x0 = x3;
		}
		fx0Val = func(x0);
		fx2Val = func(x2);
		// Print intermediate results
		if (showStep == 1)
			printf("Iteration %d, Result: %10.13f\n", k, x2);
	}
	printf("Number of iterations: %d\nResult: %10.13f\n", k, x2);
}

/* ===== ===== ===== ===== ===== ===== ===== =====
 * The new method
 * ===== ===== ===== ===== ===== ===== ===== ===== */

void newMethod() {
	printf("- New Method\n");
	double x0, x1;
	setInit(&x0, &x1, probNum);
	int k = 0;

	double fx1Val = func(x0);
	while (k < maxItr && fabs(x1 - x0) > err) {
		k++;
		// Step 1
		double fx0Val = fx1Val;
		fx1Val = func(x1);
		double denom = fx1Val - fx0Val;
		if (fabs(denom) < err) {
			break;
		}
		double xStar = x1 - (x1 - x0) / denom * fx1Val;
		// Step 2
		x0 = x1;
		denom = fx1Val - func(xStar);
		if (fabs(denom) < err)
			break;
		x1 = x1 - (x1 - xStar) / denom * fx1Val;
		// Print intermediate results
		if (showStep == 1)
			printf("Iteration %d, Result: %10.13f\n", k, x1);
	}
	printf("Number of iterations: %d\nResult: %10.13f\n", k - 1, x1);
}

/* ===== ===== ===== ===== ===== ===== ===== =====
 * The new method with bracketing
 * ===== ===== ===== ===== ===== ===== ===== ===== */

void newMethodBracketed() {
	printf("- New Method With Bracketing\n");
	double x0, x1;
	setInit(&x0, &x1, probNum);
	int k = 0;

	double fx0Val = func(x0);
	double fx1Val = func(x1);
	while (k < maxItr && fabs(fx1Val) > err && fabs(fx0Val) > err) {
		k++;
		// Step 1
		double denom = fx1Val - fx0Val;
		if (fabs(denom) < err)
			break;
		double xStar = x1 - (x1 - x0) / denom * fx1Val;
		// Step 2
		double fxStarVal = func(xStar);
		denom = fx1Val - fxStarVal;
		if (fabs(denom) < err)
			break;
		double xNew = x1 - (x1 - xStar) / denom * fx1Val;
		// Bracketing of the root
		if ((xNew - x1) * (xNew - x0) > 0) {  // Case I
			if (fx1Val * fxStarVal < 0) {
				x0 = x1;
				x1 = xStar;
			} else {
				x1 = x0;
				x0 = xStar;
			}
		} else { // Case II
			double fxNewVal = func(xNew);
			if (fxNewVal * fxStarVal < 0) {
				x0 = xNew;
				x1 = xStar;
			} else if (fxNewVal * fx0Val < 0) {
				x1 = x0;
				x0 = xNew;
			} else {
				x0 = x1;
				x1 = xNew;
			}
		}
		fx0Val = func(x0);
		fx1Val = func(x1);
		// Print intermediate results
		if (showStep)
			printf("Iteration %d, Result: %10.13f, %10.13f\n", k, x0, x1);
	}
	double result = x0;
	if (fabs(fx0Val) > fabs(fx1Val))
		result = x1;
	printf("Number of iterations: %d\nResult: %10.13f\n", k, result);
}

/* ===== ===== ===== ===== ===== ===== ===== =====
 * The secant method with Aitken's Delta 2
 * ===== ===== ===== ===== ===== ===== ===== ===== */

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
	printf("Number of iterations: %d\nResult: %10.13f\n", k, xHat);
}

/* ===== ===== ===== ===== ===== ===== ===== =====
 * Bisection method
 * ===== ===== ===== ===== ===== ===== ===== ===== */

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
	printf("Number of iterations: %d\nResult: %10.13f\n", k, (x0 + x1) * 0.5);
}

/* ===== ===== ===== ===== ===== ===== ===== =====
 * The inverse quadratic interpolation method
 * ===== ===== ===== ===== ===== ===== ===== ===== */

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
	printf("Number of iterations: %d\nResult: %10.13f\n", k, x2);
}

int main() {
	selectProb();
	initFunc();
	printf("===== ===== ===== ===== ===== ===== =====");
	printf("\nProblem %d\n\n", probNum);
	clock_t start, end;
	double cpu_time_used;

	// Newton's method
	start = clock();
	newton();
	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC * 1e6;
	printf("Time used: %5.4f microseconds\n\n", cpu_time_used);

	// Secant method
	start = clock();
	secant();
	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC * 1e6;
	printf("Time used: %5.7f microseconds\n\n", cpu_time_used);

	// SecantRF method
	start = clock();
	secantRF();
	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC * 1e6;
	printf("Time used: %5.7f microseconds\n\n", cpu_time_used);

	// Ridders' method
	start = clock();
	ridders();
	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC * 1e6;
	printf("Time used: %5.7f microseconds\n\n", cpu_time_used);

	// The new method
	start = clock();
	newMethod();
	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC * 1e6;
	printf("Time used: %5.7f microseconds\n\n", cpu_time_used);

	// The bracketed new method
	start = clock();
	newMethodBracketed();
	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC * 1e6;
	printf("Time used: %5.7f microseconds\n\n", cpu_time_used);

	// Secant method with Aitken's Delta Squared
	start = clock();
	secantAitkensDelta2();
	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC * 1e6;
	printf("Time used: %5.7f microseconds\n\n", cpu_time_used);

	// Bisection method
	start = clock();
	bisection();
	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC * 1e6;
	printf("Time used: %5.7f microseconds\n\n", cpu_time_used);

	// Inverse quadratic interpolation method
	start = clock();
	IQI();
	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC * 1e6;
	printf("Time used: %5.7f microseconds\n\n", cpu_time_used);
}
