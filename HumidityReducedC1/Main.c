#include "Fluxes.h"

int main() {
	buildMesh();
	setInitialCond();
	calcFluxes();
}
