#include "LinearMethods.h"
#include <iostream>
double dichotomySearch(double (*fncPtr)(double), double a, double b, double eps, unsigned maxIterationCount) {
	double x, y1, y2;
	//std::cout << a << " " << b << std::endl;
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
		//std::cout << a << " " << b << std::endl;
	}

	return (b + a) / 2;
}

double goldenSectionSearch(double(*fncPtr)(double), double a, double b, double eps, unsigned maxIterationCount)
{
	double goldenNumber = (pow(5, 0.5) + 1) / 2;
	double x1, x2, y1, y2;
	
	x1 = b - (b - a) / goldenNumber;
	x2 = a + (b - a) / goldenNumber;
	y1 = fncPtr(x1);
	y2 = fncPtr(x2);
	//std::cout << a << " " << b << std::endl;
	while ((b - a) / 2 > eps) {
		if (y1 < y2) {
			a = x1;
			x1 = x2;
			y1 = y2;
			x2 = a + (b - a)/goldenNumber;
			y2 = fncPtr(x2);
		}
		else {
			b = x2;
			x2 = x1;
			y2 = y1;
			x1 = b - (b - a)/goldenNumber;
			y1 = fncPtr(x1);
		}
		//std::cout << a << " " << b << std::endl;
	}
	return (a + b) / 2;
}
