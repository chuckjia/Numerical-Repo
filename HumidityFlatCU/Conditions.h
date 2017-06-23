/*
 * Conditions.h
 *
 *  Created on: Jun 21, 2017
 *      Author: chuckjia
 */

#ifndef CONDITIONS_H_
#define CONDITIONS_H_

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Initial Conditions
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

/*
 * Initialize the solution
 */
double sl[numCellsX][numCellsP][2];

/*
 * Set initial conditions
 */
void setInitCond() {
	for (int i = 0; i < numCellsX; i++)
		for (int j = 0; j < numCellsP; j++) {
			sl[i][j][0] = 0;
			sl[i][j][1] = 0;
		}
}

/* ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
 * Source Terms
 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

/*
 * Source function (component 1): wrapper
 */
double sourceFcn1(double T, double q, double x, double p, double t) {
	return 0;
}

/*
 * Source function (component 2): wrapper
 */
double sourceFcn2(double T, double q, double x, double p, double t) {
	return 0;
}

#endif /* CONDITIONS_H_ */
