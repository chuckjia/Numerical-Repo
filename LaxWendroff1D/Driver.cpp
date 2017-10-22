#include "FV.h"

int main() {
	if (methodNum == 0)
		upwind_right();
	else if (methodNum == 1)
		laxWendroff();
	else if (methodNum == 2)
		laxWendroff_limiter();
	writeResToFile();
	printMsg();
}








