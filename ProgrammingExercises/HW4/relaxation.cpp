#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <omp.h>
#define PI 3.14159265358979323846
#define TOL 1.0E-8

using namespace std;

int main()
{
	int i, N = 0;
	for(N=0; N<10000; N*=10)
	{
		double* T = new double[2*N+1];
		T[i] = 0.0;
		for(i=0; i<2*N; i++)
		{
			T[i] = 
		}
	}

	return 0;
}