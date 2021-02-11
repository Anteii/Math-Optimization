#include <iostream>
#include <cmath>
#include "LinearMethods.h"

#define PI 3.14159265358979323846
using namespace std;

double testFunc(double x) {
	return pow(x, 2);
}

int main() {
	double maxX, maxY;

	double a = 2;
	double b = 5;

	maxX = dichotomySearch(testFunc, a, b, 0.0001);
	maxY = testFunc(maxX);
	cout << maxX << " " << maxY << endl;
	maxX = goldenSectionSearch(testFunc, a, b, 0.0001);
	maxY = testFunc(maxX);
	cout << maxX << " " << maxY << endl;
	return 0;
}