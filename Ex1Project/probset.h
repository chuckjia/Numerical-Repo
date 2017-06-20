/*
 * probset.h
 *
 *  Created on: Jun 15, 2017
 *      Author: chuckjia
 */

#ifndef PROBSET_H_
#define PROBSET_H_

#include <math.h>

// Initiation of the problem number. It will be reset later.
int probNum = 1;
// The error for the stopping condition
double err = 1e-13;

/* ===== ===== ===== ===== ===== ===== ===== =====
 * Problem 1
 * ===== ===== ===== ===== ===== ===== ===== ===== */

double funcProb1(double x) {
	return pow(sin(x), 2) - x * x + 1;
}

double derFunc1(double x) {
	return 2 * sin(x) * cos(x) - 2 * x;
}

void init1(double *x0, double *x1) {
	*x0 = 1;
	*x1 = 3;
}

/* ===== ===== ===== ===== ===== ===== ===== =====
 * Problem 2
 * ===== ===== ===== ===== ===== ===== ===== ===== */

double funcProb2(double x) {
	return pow(sin(x), 2) - x * x + 1;
}

double derFunc2(double x) {
	return 2 * sin(x) * cos(x) - 2 * x;
}

void init2(double *x0, double *x1) {
	*x0 = 3;
	*x1 = 1;
}

/* ===== ===== ===== ===== ===== ===== ===== =====
 * Problem 3
 * ===== ===== ===== ===== ===== ===== ===== ===== */

double funcProb3(double x) {
	return pow(x, 2) - exp(x) - 3 * x + 2;
}

double derFunc3(double x) {
	return 2 * x - exp(x) - 3;
}

void init3(double *x0, double *x1) {
	*x0 = -50000000;
	*x1 = 3;
}

/* ===== ===== ===== ===== ===== ===== ===== =====
 * Problem 4
 * ===== ===== ===== ===== ===== ===== ===== ===== */

double funcProb4(double x) {
	return pow(x, 2) - exp(x) - 3 * x + 2;
}

double derFunc4(double x) {
	return 2 * x - exp(x) - 3;
}

void init4(double *x0, double *x1) {
	*x0 = 3;
	*x1 = -50000000;
}

/* ===== ===== ===== ===== ===== ===== ===== =====
 * Problem 5
 * ===== ===== ===== ===== ===== ===== ===== ===== */

double funcProb5(double x) {
	return x * exp(x) - 10;
}

double derFunc5(double x) {
	return (x + 1) * exp(x);
}

void init5(double *x0, double *x1) {
	*x0 = 0;
	*x1 = 2;
}

/* ===== ===== ===== ===== ===== ===== ===== =====
 * Problem 6
 * ===== ===== ===== ===== ===== ===== ===== ===== */

double funcProb6(double x) {
	return cos(x / 180 * M_PI);
}

double derFunc6(double x) {
	return - sin(x / 180 * M_PI) / 180 * M_PI;
}

void init6(double *x0, double *x1) {
	*x0 = 100;
	*x1 = 280;
}

/* ===== ===== ===== ===== ===== ===== ===== =====
 * Problem 7
 * ===== ===== ===== ===== ===== ===== ===== ===== */

double funcProb7(double x) {
	return sin(x / 180 * M_PI);
}

double derFunc7(double x) {
	return cos(x / 180 * M_PI) / 180 * M_PI;
}

void init7(double *x0, double *x1) {
	*x0 = 10;
	*x1 = 280;
}

/* ===== ===== ===== ===== ===== ===== ===== =====
 * Problem 8
 * ===== ===== ===== ===== ===== ===== ===== ===== */

double funcProb8(double x) {
	return pow(x, 3) - 2 * x - 5;
}

double derFunc8(double x) {
	return 3 * x * x - 2;
}

void init8(double *x0, double *x1) {
	*x0 = 2.5;
	*x1 = 0.01;
}

/* ===== ===== ===== ===== ===== ===== ===== =====
 * Problem 9
 * ===== ===== ===== ===== ===== ===== ===== ===== */

double funcProb9(double x) {
	return tan(x) + x - 10;
}

double derFunc9(double x) {
	return 1 / pow(cos(x), 2) + 1;
}

void init9(double *x0, double *x1) {
	*x0 = 9;
	*x1 = 10;
}

/* ===== ===== ===== ===== ===== ===== ===== =====
 * Problem 10
 * ===== ===== ===== ===== ===== ===== ===== ===== */

double funcProb10(double x) {
	return exp(tan(x)) + pow(x, 3);
}

double derFunc10(double x) {
	return 1 / pow(cos(x), 2) * exp(tan(x)) + 3 * x * x;
}

void init10(double *x0, double *x1) {
	*x0 = 1;
	*x1 = -1;
}

/* ===== ===== ===== ===== ===== ===== ===== =====
 * Problem 11
 * ===== ===== ===== ===== ===== ===== ===== ===== */

double funcProb11(double x) {
	return pow(x, 5) + pow(x, 3) + 20;
}

double derFunc11(double x) {
	return 5 * pow(x, 4) + 3 * x * x;
}

void init11(double *x0, double *x1) {
	*x0 = 2;
	*x1 = -2;
}

/* ===== ===== ===== ===== ===== ===== ===== =====
 * Problem 12
 * ===== ===== ===== ===== ===== ===== ===== ===== */

double funcProb12(double x) {
	return sin(x) * exp(x);
}

double derFunc12(double x) {
	return (cos(x) + sin(x)) * exp(x);
}

void init12(double *x0, double *x1) {
	*x0 = 2;
	*x1 = 4;
}

/* ===== ===== ===== ===== ===== ===== ===== =====
 * Problem 13
 * ===== ===== ===== ===== ===== ===== ===== ===== */

double funcProb13(double x) {
	return pow(x, 5) + pow(x, 3) + 20;
}

double derFunc13(double x) {
	return 5 * pow(x, 4) + 3 * x * x;
}

void init13(double *x0, double *x1) {
	*x0 = 2000;
	*x1 = -2000;
}

/* ===== ===== ===== ===== ===== ===== ===== =====
 * Problem 14
 * ===== ===== ===== ===== ===== ===== ===== ===== */

double funcProb14(double x) {
	return x - exp(- x * x);
}

double derFunc14(double x) {
	return 1 + (2 * x) * exp(- x * x);
}

void init14(double *x0, double *x1) {
	*x0 = 2000;
	*x1 = -2000;
}

/* ===== ===== ===== ===== ===== ===== ===== =====
 * Problem 15
 * ===== ===== ===== ===== ===== ===== ===== ===== */

double funcProb15(double x) {
	return cos(exp(x));
}

double derFunc15(double x) {
	return - exp(x) * sin(exp(x));
}

void init15(double *x0, double *x1) {
	*x0 = -1;
	*x1 = 0.5;
}

/* ===== ===== ===== ===== ===== ===== ===== =====
 * Problem 16
 * ===== ===== ===== ===== ===== ===== ===== ===== */

double funcProb16(double x) {
	return x * log(x) - 3;
}

double derFunc16(double x) {
	return log(x) + 1;
}

void init16(double *x0, double *x1) {
	*x0 = 1;
	*x1 = 1000;
}

/* ===== ===== ===== ===== ===== ===== ===== =====
 * Set Problem
 * ===== ===== ===== ===== ===== ===== ===== ===== */

double (*funcPtr) (double);
double (*derFuncPtr) (double);

void initFunc() {
	if (probNum == 1) {
		funcPtr = &funcProb1;
		derFuncPtr = &derFunc1;
	} else if (probNum == 2) {
		funcPtr = &funcProb2;
		derFuncPtr = &derFunc2;
	} else if (probNum == 3) {
		funcPtr = &funcProb3;
		derFuncPtr = &derFunc3;
	} else if (probNum == 4) {
		funcPtr = &funcProb4;
		derFuncPtr = &derFunc4;
	} else if (probNum == 5) {
		funcPtr = &funcProb5;
		derFuncPtr = &derFunc5;
	} else if (probNum == 6) {
		funcPtr = &funcProb6;
		derFuncPtr = &derFunc6;
	} else if (probNum == 7) {
		funcPtr = &funcProb7;
		derFuncPtr = &derFunc7;
	} else if (probNum == 8) {
		funcPtr = &funcProb8;
		derFuncPtr = &derFunc8;
	} else if (probNum == 9) {
		funcPtr = &funcProb9;
		derFuncPtr = &derFunc9;
	} else if (probNum == 10) {
		funcPtr = &funcProb10;
		derFuncPtr = &derFunc10;
	} else if (probNum == 11) {
		funcPtr = &funcProb11;
		derFuncPtr = &derFunc11;
	} else if (probNum == 12) {
		funcPtr = &funcProb12;
		derFuncPtr = &derFunc12;
	} else if (probNum == 13) {
		funcPtr = &funcProb13;
		derFuncPtr = &derFunc13;
	} else if (probNum == 14) {
		funcPtr = &funcProb14;
		derFuncPtr = &derFunc14;
	} else if (probNum == 15) {
		funcPtr = &funcProb15;
		derFuncPtr = &derFunc15;
	} else if (probNum == 16) {
		funcPtr = &funcProb16;
		derFuncPtr = &derFunc16;
	}
}

void setInit(double *x0, double *x1, int probNum) {
	if (probNum == 1)
		init1(x0, x1);
	else if (probNum == 2)
		init2(x0, x1);
	else if (probNum == 3)
		init3(x0, x1);
	else if (probNum == 4)
		init4(x0, x1);
	else if (probNum == 5)
		init5(x0, x1);
	else if (probNum == 6)
		init6(x0, x1);
	else if (probNum == 7)
		init7(x0, x1);
	else if (probNum == 8)
		init8(x0, x1);
	else if (probNum == 9)
		init9(x0, x1);
	else if (probNum == 10)
		init10(x0, x1);
	else if (probNum == 11)
		init11(x0, x1);
	else if (probNum == 12)
		init12(x0, x1);
	else if (probNum == 13)
		init13(x0, x1);
	else if (probNum == 14)
		init14(x0, x1);
	else if (probNum == 15)
		init15(x0, x1);
	else if (probNum == 16)
		init16(x0, x1);
	/*else if (probNum == 17)
		init17(x0, x1);
	else if (probNum == 18)
		init18(x0, x1);
	else if (probNum == 19)
		init19(x0, x1);
	else if (probNum == 20)
		init20(x0, x1);*/
}

double func(double x) {
	return (*funcPtr)(x);
}

double derFunc(double x) {
	return (*derFuncPtr)(x);
}

void selectProb() {
	printf("Enter Problem Number: ");
	scanf("%d",&probNum);
}

double sgnFunc(double x) {
	if (x > 0)
		return 1;
	else if (x < 0)
		return -1;
	return 0;
}
#endif /* PROBSET_H_ */
