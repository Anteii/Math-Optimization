#include "util.h"
#include <cmath>


double* vec_ml_vec(double* vec1, double* vec2, size_t size, double* result) {
	
	double* vec;
	if (result == NULL) {
		vec = new double[size];
	}
	else {
		vec = result;
	}
	for (size_t i = 0; i < size; ++i) {
		vec[i] = vec1[i] * vec2[i];
	}
	return vec;
}

double* vec_add_vec(double* vec1, double* vec2, size_t size, double* result){
	double* vec;
	if (result == NULL) {
		vec = new double[size];
	}
	else {
		vec = result;
	}
	for (size_t i = 0; i < size; ++i) {
		vec[i] = vec1[i] + vec2[i];
	}
	return vec;
}

void vec_add_sc(double* vec1, double a, size_t size)
{
	for (size_t i = 0; i < size; ++i) {
		vec1[i] += a;
	}
}

void vec_ml_sc(double* vec1, double a, size_t size, double* result)
{
	double* vec;
	if (result == NULL) {
		vec = vec1;
	}
	else {
		vec = result;
	}
	for (size_t i = 0; i < size; ++i) {
		vec[i] = vec1[i] * a;
	}
}

double L2norm(double* vec, size_t size) {
	double t = 0;
	for (int i = 0; i < size; ++i)
	{
		t += pow(vec[i], 2);
	}
	return pow(t, 0.5);
}

double L2distance(double* vec1, double* vec2, size_t size)
{
	double t = 0;
	for (int i = 0; i < size; ++i)
	{
		t += pow(vec1[i] - vec2[i], 2);
	}
	return pow(t, 0.5);
}

double numerical_differentiation(double(*fncPtr)(double), double x)
{
	double h = 0.0001;
	return (fncPtr(x + h) - fncPtr(x - h)) / h / 2;
}

double numerical_differentiation(double(*fncPtr)(double*, size_t n), double* v, size_t k, size_t n)
{
	double y1, y2;
	double h = 0.00001;
	double* temp = new double[n];
	for (int i = 0; i < n; ++i) temp[i] = v[i];
	temp[k] -= h;
	y1 = fncPtr(temp, n);
	temp[k] += 2 * h;
	y2 = fncPtr(temp, n);
	delete[] temp;
	return (y2 - y1) / (2 * h);
}

double* numerical_gradient(double (*fncPtr)(double*, size_t n), double* v, size_t n) {
	double* grad = new double[n];
	for (size_t i = 0; i < n; ++i) {
		grad[i] = numerical_differentiation(fncPtr, v, i, n);
	}
	return grad;
}

double calc_mixed_derivative(double(*fncPtr)(double*, size_t n), double* v, size_t k1, size_t k2, size_t n) {
	double* nv = new double[n];
	double* y = new double[4];
	double h1 = 1e-5, h2 = 1e-5;

	for (int i = 0; i < n; i++) nv[i] = v[i];

	nv[k1] += h1;
	nv[k2] += h2;
	y[0] = fncPtr(nv, n);

	nv[k2] -= 2 * h2;
	y[1] = fncPtr(nv, n);

	nv[k1] -= 2*h1;
	nv[k2] += 2*h2;
	y[2] = fncPtr(nv, n);

	nv[k2] -= 2 * h2;
	y[3] = fncPtr(nv, n);
	return (y[0] - y[1] - y[2] + y[3]) / (4 * h1 * h2);
}

double** calc_hessian(double(*fncPtr)(double*, size_t n), double* v, size_t k, size_t n)
{
	double** hessian = new double* [n];
	for (int i = 0; i < n; i++)
	{
		hessian[i] = new double[n];
		for (int j = i; i < n; i++)
		{
			hessian[i][j] = calc_mixed_derivative(fncPtr, v, i, i, n);
			hessian[j][i] = hessian[i][j];
		}

	}
	return hessian;
}

