/*
 * InitialCondtions.h
 *
 *  Created on: Jun 19, 2017
 *      Author: chuckjia
 */

#ifndef INITIALCONDITIONS_H_
#define INITIALCONDITIONS_H_

#include "Mesh.h"

void setInitialCond() {
	for (int i = 0; i < numCellsXDir; i++)
		for (int j = 0; j < numCellsPDir; j++) {
			soln[i][j][0] = 1;
			soln[i][j][1] = 0;
		}
}


#endif /* INITIALCONDITIONS_H_ */
