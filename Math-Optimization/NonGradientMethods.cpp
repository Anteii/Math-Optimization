#include "NonGradientMethods.h"
#include "util.h"
#include "LinearMethods.h"


double* gauss_zeidel_optimization(double(*fncPtr)(double*, size_t), double* x0, size_t n, double eps)
{
	double* x_prev = new double[n];
	double* x_cur = new double[n];
	double* lambdas = new double[n];
	double delta;

	for (size_t i = 0; i < n; ++i) {
		lambdas[i] = 0;
	}
	for (size_t i = 0; i < n; ++i) {
		x_prev[i] = x0[i];
		x_cur[i] = x0[i];
	}

	int i = 0;
	do {
		x_prev[i] = x_cur[i];
		delta = x_prev[i]/2 + 1;

		lambdas[i] = dichotomy_min_search(fncPtr, x_prev, i, n, -delta, delta, eps);
		
		x_cur[i] += lambdas[i];

		// counter up
		++i;
		if (i == n) i = 0;
	} while (L2distance(x_cur, x_prev, n) > eps);

	delete[] x_prev;
	delete[] lambdas;

	return x_cur;
}
