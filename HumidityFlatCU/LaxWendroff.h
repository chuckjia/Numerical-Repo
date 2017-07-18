/*
 * LaxWendroff.h
 *
 *  Created on: Jul 5, 2017
 *      Author: chuckjia
 */

#ifndef LAXWENDROFF_H_
#define LAXWENDROFF_H_

double DtSq = Dt * Dt;
double DxInvSq = DxInv * DxInv;
double DpInvSq = DpInv * DpInv;
double DpSq = Dp * Dp;

void lwOneStep() {
	for (int i = 1; i < lastIndexX; i++)
		for (int j = 1; j < lastIndexP; j++)
			for (int ii = 0; ii < 2; ii++) {
				double A = getuFcnVal(i, j), B = getomegaFcnVal(i, j);
				sl[i][j][ii] += 0.5 *
						(- Dt * DxInv * A * (sl[i+1][j][ii] - sl[i - 1][j][ii])
								- Dt * DpInv * B * (sl[i][j+1][ii] - sl[i][j - 1][ii])
								+ DtSq * DxInvSq * A * A * (sl[i + 1][j][ii] - 2 * sl[i][j][ii] + sl[i - 1][j][ii])
								+ DtSq * DpInvSq * B * B * (sl[i][j + 1][ii] - 2 * sl[i][j][ii] + sl[i][j - 1][ii])
								+ 0.25 * DtSq * DxInv * DpInv * 2 * A * B *
								(sl[i + 1][j + 1][ii] - sl[i - 1][j + 1][ii] - (sl[i + 1][j - 1][ii] - sl[i - 1][j - 1][ii])));
			}
}

void lwSourceStep(int i, int j) {

}

void lw() {
	printf("Using Lax-Wendroff Method");
	for (int tt = 0; tt < numTimeSteps; tt++) {
		// Temporary value for the progress message
		double timer_factor1_CONST = 100.0 / finalTime;
		double t = tt * Dt;
		// Printing progress
		printf("\r  [Progress]:%5.1f%%", t * timer_factor1_CONST);
		fflush(stdout);
		if (tt % 2 == 0)
			lwOneStep();
		else {
			for (int i = 1; i < lastIndexX; i++)
				for (int j = 1; j < lastIndexP; j++) {
					double T = sl[i][j][0], q = sl[i][j][1];
					double sourceVal[2] = {0, 0};
					double x = getCellCenterX(i, j), p = getCellCenterP(i, j);
					(*sourceFcnPtr)(sourceVal, T, q, x, p, t, i, j);
					double S1[2] = {0, 0};
					S1[0] = T + 0.5 * Dt * sourceVal[0];
					S1[1] = q + 0.5 * Dt * sourceVal[0];

					sourceVal[0] = 0; sourceVal[1] = 0;
					(*sourceFcnPtr)(sourceVal, S1[0], S1[1], x, p, t + 0.5 * Dt, i, j);
					for (int jj = 0; jj < 2; jj++)
						sl[i][j][jj] += Dt * sourceVal[jj];
				}
		}
		dirichletCond();
	}
	printf("\r [Complete]: 100%%   \n");
}

#endif /* LAXWENDROFF_H_ */
