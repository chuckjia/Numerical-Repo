/*
 * Driver.cpp
 *
 *  Created on: Mar 13, 2018
 *      Author: chuckjia
 */

#include <stdio.h>
#include <math.h>
#include "../include/ModelParam.h"
#include "MeshBuild.h"
#include "../include/debug.h"
#include "../include/NumSoln.h"
#include "../Models/Model_0.h"
#include "Testing.h"

void print_msg(int k) {
	if (k == 0) {
		// Starting message
		printf("===== ===== ===== ===== ===== ===== ===== ===== ===== ===== \n");
		printf(">> Numerical Simulation of the Primitive Equation Model\n");
		printf("===== ===== ===== ===== ===== ===== ===== ===== ===== ===== \n\n");
	} else if (k == 1) {
		// Ending message
		printf("===== ===== ===== ===== ===== ===== ===== ===== ===== ===== \n");
		printf(">> Simulation Completed\n");
		printf("===== ===== ===== ===== ===== ===== ===== ===== ===== ===== \n\n");
	}
}

void numer_simulation(ModelParam &model) {

	print_msg(0);

	/*
	 * Build mesh
	 */

	int num_cell_x = model.numCellX(), num_cell_p = model.numCellP(), num_gridpt_x = num_cell_x + 1, num_gridpt_p = num_cell_p + 1;

	// Build mesh
	MeshGrid grid_x(num_gridpt_x, num_gridpt_p), grid_p(num_gridpt_x, num_gridpt_p);
	Matrix2D cell_centers_x(num_cell_x, num_cell_p), cell_centers_p(num_cell_x, num_cell_p);
	Vector1D cell_left_Dp(num_gridpt_x);
	MeshBuild::buildMesh(model, grid_x, grid_p, cell_centers_x, cell_centers_p, cell_left_Dp);

	/*
	 * Enforce initial conditions
	 */

	NumSoln T(num_cell_x, num_cell_p), q(num_cell_x, num_cell_p), u(num_cell_x, num_cell_p),
			w(num_cell_x, num_cell_p), phi_x(num_cell_x, num_cell_p);


	print_msg(1);
}

int main() {
	// ModelParam model = create_model_0(10, 10);
	// numer_simulation(model);
	test();
}
