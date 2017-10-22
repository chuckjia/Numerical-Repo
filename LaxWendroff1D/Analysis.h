/*
 * Analysis.h
 *
 *  Created on: Oct 8, 2017
 *      Author: chuckjia
 */

#ifndef ANALYSIS_H_
#define ANALYSIS_H_
#include "Models.h"

void writeResToFile() {
	FILE *f = fopen("res.txt", "wb");
	for (int i = 1; i < lastGhostIndex; i++)
		fprintf(f, "%f ", u[i]);
	fclose(f);
	FILE *g = fopen("par.txt", "wb");
	fprintf(g, "%f ", Dt);
	fprintf(g, "%d", numTimeSteps);
	fclose(g);
}

void writeResToFileForMovie(int tt) {
	char filename[20];
	sprintf(filename, "res%d.txt", tt);
	FILE *f = fopen(filename, "wb");
	for (int i = 1; i < lastGhostIndex; i++)
		fprintf(f, "%f ", u[i]);
	fclose(f);
}

double calcL2Err() {
	double l2Err = 0;
	for (int i = 1; i < lastGhostIndex; i++) {
		double x = getCellCenter(i);
		double exactSolnVal = (*exactSolnPtr)(x, finalTime);
		l2Err += pow(u[i] - exactSolnVal, 2);
	}
	l2Err = sqrt(Dx * l2Err);
	return l2Err;
}

void printMsg() {
	printf("\n");
	if (methodNum == 0)
		printf("Applying Upwind-type Godunov Method\n");
	else if (methodNum == 1)
		printf("Applying Lax-Wendroff Method without Flux Limiters\n");
	else if (methodNum == 2)
			printf("Applying Lax-Wendroff Method with Flux Limiters\n");
	printf("  Model %d\n", modelNum);
	printf("  - Dt = %1.4es, finalTime = %1.2es\n", Dt, finalTime);
	printf("  - Nx = %d, Dx = %1.4e, Dt/Dx = %1.2f\n", numDivisions, Dx, Dt / Dx);
	printf("  - L2 Err = %1.7e", calcL2Err());
	printf("\n");
}

#endif /* ANALYSIS_H_ */
