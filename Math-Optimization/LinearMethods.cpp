#include "LinearMethods.h"
#include "util.h"

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
	Dichotomy method for function f(X), where x = (x1, x2, .., lambda + xk, .., xn)
	Minimize along k axis
	fncPtr - function pointer
	v - value vector
	k - optimizing argument position
	n - dimension number
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
/*
	Dichotomy method for function f(X), where X = v + lambda * w
	fncPtr - function pointer
	v - value vector
	w - weights
	n - dimension number
	a - left interval point
	b - right interval point
	eps - precision
*/
double dichotomy_min_search(double (*fncPtr)(double*, size_t n), double* v, double* w, size_t n, double a, double b, double eps) {

	double x, y1, y2;
	double* term_1 = new double[n];
	double* grad = new double[n];

	while (b - a > eps) {
		x = (b + a) / 2 - eps / 2;

		vec_ml_sc(w, x, n, grad);
		vec_add_vec(grad, v, n, term_1);
		y1 = fncPtr(term_1, n);

		x += eps;
		vec_ml_sc(w, x, n, grad);
		term_1 = vec_add_vec(grad, v, n, term_1);
		y2 = fncPtr(term_1, n);
		
		x -= eps / 2;
		
		if (y1 < y2) {
			b = x;
		}
		else {
			a = x;
		}
	}
	delete[] grad;
	delete[] term_1;
	return (b + a) / 2;
}