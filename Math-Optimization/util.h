#pragma once
#include <string>
double* vec_add_vec(double* vec1, double* vec2, size_t size, double* result = NULL);
void vec_add_sc(double* vec1, double a, size_t size);
void vec_ml_sc(double* vec1, double a, size_t size, double* result = NULL);
double* vec_ml_vec(double* vec1, double* vec2, size_t size, double* result = NULL);

double L2norm(double* vec, size_t size);
double L2distance(double* vec1, double* vec2, size_t size);

double numerical_differentiation(double (*fncPtr)(double), double x);
double numerical_differentiation(double (*fncPtr)(double*, size_t n), double* v, size_t k, size_t n);
double* numerical_gradient(double (*fncPtr)(double*, size_t n), double* v, size_t n);