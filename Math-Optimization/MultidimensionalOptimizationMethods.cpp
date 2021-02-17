#include "MultidimensionalOptimizationMethods.h"
#include "util.h"
#include "LinearMethods.h"


double* gauss_zeidel_optimization(double(*fncPtr)(double*, size_t), double* x0, size_t n, double eps)
{
	double* x_prev = new double[n];
	double* x_cur = new double[n];
	double* lambdas = new double[n];
	double delta;

	for (size_t j = 0; j < n; ++j) {
		x_prev[j] = x0[j];
		x_cur[j] = x0[j];
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
	} while (L2norm(x_cur, x_prev, n) > eps);

	delete[] x_prev;
	delete[] lambdas;

	return x_cur;
}
