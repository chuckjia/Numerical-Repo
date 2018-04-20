/*
 * Model_0.h
 *
 *  Created on: Mar 20, 2018
 *      Author: chuckjia
 */

#ifndef SRC_MODEL_0_H_
#define SRC_MODEL_0_H_

#include <math.h>
#include "../include/ModelParam.h"

double pB_MDL0(double x) {
	double term = (x - 37500) / 6000.;
	return 1000. - 250. * exp(-term * term);
}

double pB_x_MDL0(double x) {
	double term = (x - 37500) / 6000.;
	return 1 / 12. * term * exp(-term * term);
}

double init_T_MDL0(double x, double p, double t) {
	return 300 - (1 - p * 0.001) * 50;
}

double es_fcn_MDL0(double T) {
	return 6.112 * exp(17.67 * (T - 273.15) / (T - 29.65));
}

double init_q_MDL0(double x, double p, double t) {
	double T = 300 - (1 - p * 0.001) * 50;
	return 0.622 * es_fcn_MDL0(T) / p;
}

double init_u_MDL0(double x, double p, double t) {
	int n = 1;
	return 7.5 + 2 * cos(p * M_PI * 0.001) * cos(2 * n * M_PI * x / 75000);
}

ModelParam create_model_0() {
	int Nx = 100, Np = 100;
	ModelParam model(0, 75000, 200, pB_MDL0, pB_x_MDL0, Nx, Np, init_T_MDL0, init_q_MDL0, init_u_MDL0);
	return model;
}

ModelParam create_model_0(int Nx, int Np) {
	ModelParam model(0, 75000, 200, pB_MDL0, pB_x_MDL0, Nx, Np, init_T_MDL0, init_q_MDL0, init_u_MDL0);
	return model;
}


#endif /* SRC_MODEL_0_H_ */
