#pragma once
#include "GradientMethods.h"
#include "LinearMethods.h"
#include "util.h"


double* gradient_method_optimization(double(*fncPtr)(double*, size_t), double* x0, size_t n, double eps)
{
	double lambda;
	double delta;
	double* x_cur = new double[n];
	double* x_prev = new double[n];
	double* grad;

	for (size_t i = 0; i < n; ++i) {
		x_cur[i] = x0[i];
	}

	do {
		// x_prev = x_cur
		memcpy(x_prev, x_cur, sizeof(double) * n);
		
		delta = L2norm(x_prev, n)/2 + 1;
		grad = numerical_gradient(fncPtr, x_prev, n);
		
		lambda = dichotomy_min_search(fncPtr, x_prev, grad, n, -delta, delta, eps);
		// Calculate x_cur
		vec_ml_sc(grad, lambda, n);
		delete[] x_cur;
		x_cur = vec_add_vec(grad, x_prev, n);
		delete[] grad;
		double tt = fncPtr(x_prev, n);
	} while( L2distance(x_cur, x_prev, n) > eps);
	delete[] x_prev;
	return x_cur;
}
