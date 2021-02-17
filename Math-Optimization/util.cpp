#include "util.h"
#include <cmath>

double L2norm(double* vec1, double* vec2, size_t size)
{
	double t = 0;
	for (int i = 0; i < size; ++i)
	{
		t += pow(vec1[i] - vec2[i], 2);
	}
	return pow(t, 0.5);
}
