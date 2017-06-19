#include <stdio.h>

void testFunc(int m, int n, double arr[m][n]) {
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			arr[i][j] = 1;
}

int main() {
	double x[10][10];
	for(int i = 0; i < 10; i++)
		for (int j = 0; j < 10; j++) {
			x[i][j] = 0;
			printf("x[%d][%d] == %10.20f\n", i, j, x[i][j]);
		}

	testFunc(10, 10, x);
	for(int i = 0; i < 10; i++)
		for (int j = 0; j < 10; j++) {
			printf("x[%d][%d] == %10.20f\n", i, j, x[i][j]);
		}
}
