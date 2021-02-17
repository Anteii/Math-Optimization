#pragma once

double* gauss_zeidel_optimization(double(*fncPtr)(double*, size_t), double* x0, size_t n, double eps);
