#pragma once
#include <cmath>


double dichotomy_max_search(double (*fncPtr)(double), double a, double b, double eps);
double golden_section_max_search(double (*fncPtr)(double), double a, double b, double eps);
double fibonacci_max_search(double (*fncPtr)(double), double a, double b, double eps);
double dichotomy_min_search(double (*fncPtr)(double*, size_t n), double* v, size_t k, size_t n, double a, double b, double eps);
double dichotomy_min_search(double (*fncPtr)(double*, size_t n), double* v, double* w, size_t n, double a, double b, double eps);