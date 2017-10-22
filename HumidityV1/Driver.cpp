/*
 * Driver.cpp
 *
 *  Created on: Oct 13, 2017
 *      Author: chuckjia
 */

#include "Testing.h"

void testing() {
	// printMeshToFile();
	// printPar();
	// printDiagnostics();
	 test_GaussElimProj();
}

int main() {
	setMesh();
	testing();
}
