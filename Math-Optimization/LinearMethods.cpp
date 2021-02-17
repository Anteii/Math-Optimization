#include "LinearMethods.h"


double dichotomy_max_search(double (*fncPtr)(double), double a, double b, double eps) {
	double x, y1, y2;
	
	while (b - a > eps) {
		x = (b + a) / 2;
		y1 = fncPtr(x - eps/2);
		y2 = fncPtr(x + eps/2);
		
		if (y1 > y2) {
			b = x;
		}
		else {
			a = x;
		}
	}
	return (b + a) / 2;
}

double golden_section_max_search(double(*fncPtr)(double), double a, double b, double eps){
	double goldenNumber = (pow(5, 0.5) + 1) / 2;
	double x1, x2, y1, y2;
	
	x1 = b - (b - a) / goldenNumber;
	x2 = a + (b - a) / goldenNumber;
	y1 = fncPtr(x1);
	y2 = fncPtr(x2);

	while (b - a > eps) {
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

	}
	return (a + b) / 2;
}

double fibonacci_max_search(double(*fncPtr)(double), double a, double b, double eps){
	double fib[3] = { 1, 2, 3 };
	double x1, x2;
	double Ln = (b - a) / eps;

	while (fib[2] < Ln) {
		x1 = a + fib[0] / fib[2] * (b - a);
		x2 = a + fib[1] / fib[2] * (b - a);

		if (fncPtr(x1) < fncPtr(x2)) {
			a = x1;
		}
		else if (fncPtr(x1) > fncPtr(x2)) {
			b = x2;
		}

		fib[0] = fib[1];
		fib[1] = fib[2];
		fib[2] = fib[0] + fib[1];
	}
	return (a + b) / 2;
}

/*
	fncPtr - function pointer
	v - value vector
	k - optimizing argument position
	n - v dimension
	a - left interval point
	b - right interval point
	eps - precision
*/
double dichotomy_min_search(double (*fncPtr)(double*, size_t n), double* v, size_t k, size_t n, double a, double b, double eps) {

	double x, y1, y2;

	while (b - a > eps) {
		x = (b + a) / 2;
		v[k] += x;
		v[k] -= eps / 2;
		y1 = fncPtr(v, n);
		v[k] += eps;
		y2 = fncPtr(v, n);
		v[k] -= eps / 2 + x;
		if (y1 < y2) {
			b = x;
		}
		else {
			a = x;
		}
	}
	return (b + a) / 2;
}