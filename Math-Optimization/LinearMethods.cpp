#include "LinearMethods.h"
double dichotomySearch(double (*fncPtr)(double), double a, double b, double eps, unsigned maxIterationCount) {
	double x, y1, y2;
	
	while (b - a > eps) {
		x = (b + a) / 2;
		y1 = fncPtr(x - eps);
		y2 = fncPtr(x + eps);
		if (y1 > y2) {
			b = x;
		}
		else {
			a = x;
		}
	}

	return (b + a) / 2;
}

double goldenSectionSearch(double(*fncPtr)(double), double a, double b, double eps, unsigned maxIterationCount)
{
	double goldenNumber = (pow(5, 0.5) + 1) / 2;
	double x1, x2;
	
	x1 = a + (b - a) / goldenNumber;
	x2 = b - (b - a) / goldenNumber;
	while ((b - a) / 2 > eps) {
		if (fncPtr(x1) > fncPtr(x2)) {
			a = x1;
			x1 = x2;
			x2 = b - (b - a)/goldenNumber;
		}
		else {
			b = x2;
			x2 = x1;
			x1 = a + (b - a)/goldenNumber;
		}
	}
	return (a + b) / 2;
}
