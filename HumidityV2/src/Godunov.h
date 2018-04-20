/*
 * Godunov.h
 *
 *  Created on: Mar 28, 2018
 *      Author: chuckjia
 */

#ifndef SRC_GODUNOV_H_
#define SRC_GODUNOV_H_

#include "../include/MeshGrid.h"
#include "../include/ModelParam.h"

class Godunov {
public:
	void upwind(MeshGrid meshgrid, Matrix2D GG_T, Matrix2D GG_q, Matrix2D GG_u, Matrix2D FF_T, Matrix2D FF_q, Matrix2D FF_u) {

	}
};

#endif /* SRC_GODUNOV_H_ */
