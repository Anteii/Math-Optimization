#include <iostream>
#include <cmath>
#include <iomanip>
#include "util.h"
#include "LinearMethods.h"
#include "NonGradientMethods.h"
#include "GradientMethods.h"
#include "SimplexMethod.h"
#define PI 3.14159265358979323846

using namespace std;

// Function to inspect
double test_func_1(double x) {
	//return sin(pow(x, 2)) + pow(x, -45) + x;
	return x * x;
}

double test_func_2(double* x, size_t n) {
	//return x[0]*x[0]*4 + x[1]*x[1]*2;
	//return x[0] * x[0] + 1.7 * x[1] * x[1] + 2 * x[0] - 1.5 * x[1] + 0.1 * x[0] * x[1];
	return sin(x[0]) + x[1] * x[1];
}
void test_linear_methods() {
	double maxX, maxY;

	double a = 2;
	double b = 6;
	double eps = 1e-9;

	// to beautify output
	cout << setprecision(7);
	cout << "                          Xmax     Ymax" << endl;

	cout << "Dichotomy search       ";
	maxX = dichotomy_max_search(test_func_1, a, b, eps);
	maxY = test_func_1(maxX);
	cout << maxX << " " << maxY << endl;

	cout << "Golden section search  ";
	maxX = golden_section_max_search(test_func_1, a, b, eps);
	maxY = test_func_1(maxX);
	cout << maxX << " " << maxY << endl;

	cout << "Fibonacci search       ";
	maxX = fibonacci_max_search(test_func_1, a, b, eps);
	maxY = test_func_1(maxX);
	cout << maxX << " " << maxY << endl;
}

void test_gauss_zeidel() {
	size_t n = 2;
	double vec0[2] = { 0.5, 5 };
	double eps = 0.000001;

	double* x_vec = gauss_zeidel_optimization(test_func_2, vec0, n, eps);
	double min = test_func_2(x_vec, 2);

	cout << "MD optimization" << endl;
	cout << "X(";
	for (size_t i = 0; i < n; ++i) {
		cout << x_vec[i] << (i == n - 1 ? "" : "; ");
	}
	cout << ") f(X)=" << min << endl;
	delete[] x_vec;
}

void test_gradient_method() {
	size_t n = 2;
	double vec0[2] = { 0.5, 5 };
	double eps = 0.000001;

	double* x_vec = gradient_method_optimization(test_func_2, vec0, n, eps);
	double min = test_func_2(x_vec, 2);

	cout << "Gradient method optimization" << endl;
	cout << "X(";
	for (size_t i = 0; i < n; ++i) {
		cout << x_vec[i] << (i == n - 1 ? "" : "; ");
	}
	cout << ") f(X)=" << min << endl;
	delete[] x_vec;
}
double ta(double* args, size_t n) {
	return sin(args[0]) * sin(args[1]);
}

double taa(double* args, size_t n) {
	return cos(args[0]) * cos(args[1]);
}

void testSimplex() {
	std::string f = "+2x_0+1x_1";
	std::string a = "+3x_0+1x_1>=3";
	std::string b = "+4x_0+3x_1>=6";
	std::string c = "+1x_0+2x_1<=3";

	vector<double> fvec = parseExpression(f);
	Expression l1 = parseLimitation(a);
	Expression l2 = parseLimitation(b);
	Expression l3 = parseLimitation(c);

	vector<Expression> limitations({ l1, l2, l3 });
	doSimplex(fvec, limitations);
	system("pause");
}
int main() {
	//test_linear_methods();
	//test_gauss_zeidel();
	//test_gradient_method();.
	testSimplex();
	return 0;
}