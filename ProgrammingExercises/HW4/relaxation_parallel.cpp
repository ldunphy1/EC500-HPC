#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <omp.h>
#define MAX 1000000

using namespace std;

int jacobi(int N)
{
	int i = 0;
	int j = 0;
	int k = 0;
	double* T = new double[2*N+1];
	double* b = new double[2*N+1];
	double* r = new double[2*N+1];
	double* temp = new double[2*N+1];
	double num = 0.0;
	double den = 0.0;
	double normResidual = 0.0;
	for(i=0; i<2*N+1; i++)
	{
		T[i] = 0.0;
		b[i] = 0.0;
		temp[i] = 0.0;
	}
	b[N] = 1.0;
	for(k=0;k<MAX;k++)
	{
		#pragma omp parallel for shared(N, temp, T, b)
		for(i=1;i<2*N;i++)
			temp[i] = 0.5*(T[i+1]+T[i-1]) + b[i];
		#pragma omp parallel for shared(N, temp, T, b)
		for(j=1;j<2*N;j++)
			T[j] = temp[j];		
		num = 0.0;
		den = 0.0;	
		#pragma omp parallel for shared(N, T, b, r) reduction(+:num, den)	
		for(i=1;i<2*N;i++)
		{
			r[i] = b[i] - (T[i] - 0.5*(T[i-1]+T[i+1]));
			num += r[i] * r[i];
			den += b[i] * b[i];
		}
		normResidual = sqrt(num)/sqrt(den);
		if(normResidual < pow(10.0,-6.0))
		{
			printf("%f\n", normResidual);
			return k;
		}
	}
	printf("Reached MAX interations\n");
	return 0;
}

int main()
{
	int iterations = 0;
	for(int N =10; N<10000; N+=10)
	{
		iterations = jacobi(N);
		printf("%d %d\n", N, iterations);
	}
}