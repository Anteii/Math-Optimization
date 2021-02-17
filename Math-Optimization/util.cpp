#include "util.h"
#include <cmath>

double L2norm(double* vec1, double* vec2, size_t size)
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
	double h = 0.0001;
	double* temp = new double[n];
	for (int i = 0; i < n; ++i) temp[i] = v[i];
	temp[k] -= h;
	y1 = fncPtr(temp, n);
	temp[k] += 2 * h;
	y2 = fncPtr(temp, h);
	delete[] temp;
	return (y2 - y1) / h / 2;
}
