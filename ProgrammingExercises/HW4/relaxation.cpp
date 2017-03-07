#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#define MAX 100000

using namespace std;

int main()
{
	int i = 0;
	int j = 0;
	int k = 0;
	int N = 0;
	for(N=10; N<10000; N*=10)
	{
		double* T = new double[2*N+1];
		double* b = new double[2*N+1];
		double* r = new double[2*N+1];
		double* temp = new double[2*N+1];
		double num = 0.0;
		double den = 0.0;
		for(i=0; i<2*N+1; i++)
		{
			T[i] = 0.0;
			b[i] = 0.0;
			temp[i] = 0.0;
		}
		b[N] = 1.0;
		for(k=0;k<MAX;k++)
		{
			for(i=1;i<2*N;i++)
			{
				temp[i] = 0.5*(T[i+1]+T[i-1]) + b[i];
			}
			for(j=i;j<2*N;j++)
			{
				T[j] = temp[j];				
			}
			for(i=0;i<=2*N;i++)
			{
				r[i] = b[i] - (T[i] - 0.5*(T[i-1]+T[i+1]));
				num += r[i] * r[i];
				den += b[i] * b[i];
				printf("%f ", r[i]);
			}
			if((sqrt(num)/sqrt(den)) < pow(10.0,-6.0))
				return 0;
			printf("\n");
		}
	}

	return 0;
}