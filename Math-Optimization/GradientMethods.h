#pragma once

double* gradient_method_optimization(double(*fncPtr)(double*, size_t), double* x0, size_t n, double eps);
double* conjugate_gradient_method_optimization();