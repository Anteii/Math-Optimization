#pragma once
#include <cmath>
double dichotomySearch(double (*fncPtr)(double), double a, double b, double eps);
double goldenSectionSearch(double (*fncPtr)(double), double a, double b, double eps);
double fibonacciSearch(double (*fncPtr)(double), double a, double b, double eps);