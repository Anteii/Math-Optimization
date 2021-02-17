#pragma once

double L2norm(double* vec1, double* vec2, size_t size);
double numerical_differentiation(double (*fncPtr)(double), double x);
double numerical_differentiation(double (*fncPtr)(double*, size_t n), double* v, size_t k, size_t n);