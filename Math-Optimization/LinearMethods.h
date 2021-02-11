#pragma once
#include <cmath>
double dichotomySearch(double (*fncPtr)(double), double a, double b, double eps, unsigned maxIterationCount = 10000);
double goldenSectionSearch(double (*fncPtr)(double), double a, double b, double eps, unsigned maxIterationCount = 10000);