#include <iostream>
#include <cmath>
#include <iomanip>
#include "LinearMethods.h"

#define PI 3.14159265358979323846

using namespace std;

// Function to inspect
double testFunc(double x) {
	//return sin(pow(x, 2)) + pow(x, -45) + x;
	return x * x;
}

int main() {
	double maxX, maxY;
	/*
	a - left point
	b - right point
	eps - precision
	*/
	double a = 2;
	double b = 6;
	double eps = 1e-9;

	// to beautify output
	cout << setprecision(7);
	cout << "                          Xmax     Ymax" << endl;
	
	cout << "Dichotomy search       ";
	maxX = dichotomySearch(testFunc, a, b, eps);
	maxY = testFunc(maxX);
	cout << maxX << " " << maxY << endl;

	cout << "Golden section search  ";
	maxX = goldenSectionSearch(testFunc, a, b, eps);
	maxY = testFunc(maxX);
	cout << maxX << " " << maxY << endl;

	cout << "Fibonacci search       ";
	maxX = fibonacciSearch(testFunc, a, b, eps);
	maxY = testFunc(maxX);
	cout << maxX << " " << maxY << endl;

	return 0;
}