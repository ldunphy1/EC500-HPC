#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <omp.h>
#define PI 3.14159265358979323846
#define TOL 1.0E-8

using namespace std;

double f(double p0, double p1, double p2, double p3, double m)
{
	return (1/pow((pow(sin(PI*p0/2.0), 2.0) + (pow(sin(PI*p1/2.0), 2.0)) + pow(sin(PI*p2/2.0), 2.0) + pow(sin(PI*p3/2.0), 2.0) + pow(m, 2.0)), 2.0));
}

void gaussianElimination(double** A, double* b, int dim)
{ 
	int max, i, j, k = 0;
	double* tempA = new double[dim];
	for(i=0;i<dim; i++)
	{
		tempA[i] = 0.0;
	}
	double tempB, factor, sum = 0.0;
	for(i=0;i<dim;i++)
	{
		max = i;
		for(j=i+1;j<dim;j++)
		{
			if(abs(A[j][i]) > abs(A[max][i]))
			{
				max = j;
			}
		}
		tempA = A[i];
		A[i] = A[max];
		A[max] = tempA;
		tempB = b[i];
		b[i] = b[max];
		b[max] = tempB;
		for(j=i+1;j<dim;j++)
		{
			factor = A[j][i]/A[i][i];
			b[j] -= factor * b[i];
			for(k=i;k<dim;k++)
			{
				A[j][k] -= factor * A[i][k];
			}
		}
	}
	double* temp = new double[dim];
	for(i=0;i<dim;i++)
	{
		temp[i] = 0.0;
	}
	for(i=dim-1;i>=0;i--)
	{
		sum = 0.0;
		for(j=i+1;j<dim;j++)
		{
			sum += A[i][j] * temp[j];
		}
		temp[i] = (b[i] - sum) / A[i][i];
		b[i] = temp[i];
	}
	delete[] tempA, temp;
}

void getLegendreCoeff(double* a, int n)
{		
	a[0+0*(n+1)] = 1.0; //set P0
	a[1+1*(n+1)] = 1.0; //set P1
    for(int pn=2; pn<n+1;pn++)
    { 
    	a[0 +pn*(n+1)] = -((pn-1)*a[0 + (pn-2)*(n+1)])/(double) pn;
        for(int i =1;i < pn+1;i++)
    		a[i+pn*(n+1)] = ((2*pn-1)*a[i-1 +(pn-1)*(n+1)] - (pn-1)*a[i + (pn-2)*(n+1)])/(double) pn;
    }			
}

double func(double* a, double x, int n)
{
	double sum = a[0];
	for(int i=1;i<=n;i++)
	{
		sum += a[i] * pow(x,i);
	}
	return sum;
}

double fprime(double* a, double x, int n)
{
	double sum = a[1];
	for(int i=2;i<=n;i++)
	{
		sum += a[i] * i * pow(x,i-1.0);
	}
	return sum;
}

void getLegendreZero(double* zero, double* a, int n)
{
	int k, count=0;
	double xprev, xnext, f, fp = 0.0;
	bool initial=true;
	for (k=1;k<=n;k++)
	{
		initial = true;
		xprev = (1.0 - (1.0/(8.0 * pow(n,2.0))) + (1.0/(8.0 * pow(n,3.0)))) * cos(PI * (4.0*k - 1.0)/(4.0*n + 2.0));
		xnext = 2.0;
		while ( (abs((xnext-xprev)/xprev) >= TOL) && (xprev!=0.0) )
		{
			if(initial==false)
			{
				xprev = xnext;
			}
			f = func(a, xprev, n);
			fp = fprime(a, xprev, n);
			xnext = xprev - f/fp;
			initial=false;
			count++;
		}
		zero[k-1] = xnext;
	} 
}

double gaussQuad(int n, double m)
{
	int i, j, k, l = 0;
	double gqIntegral = 0.0;
	double* w = new double[n];
	double** system = new double*[n];
	double* zero = new double[n+1];
	double* a = new double[n+1];
	double* A = new double[(n+1)*(n+1)];
	for(i=0; i<n+1;i++)
	{
		zero[i] = 0.0;
		a[i] = 0.0;
	}
	for(i=0;i<(n+1)*(n+1);i++)
	{
		A[i] =0.0;
	}
	for (i=0;i<n;i++) { system[i] = new double[n]; }
	getLegendreCoeff(A, n);
	for (i=0;i<=n;i++) 
	{ 
		a[i] = A[i + n*(n+1)]; 
	}
	for(i=0;i<n;i++)
	{
		w[i] = 0.0;
		for(k=0;k<n;k++)
		{
			system[i][k] = 0.0;
		}
	}
	getLegendreZero(zero, a, n);
	for(k=0;k<n;k++)
	{
		if(k%2)
		{
			w[k] = 0.0;
		}
		else
		{
			w[k] = 2.0/(k+1.0);
		}
	}
	for(i=0;i<n;i++)
	{
		for(k=0;k<n;k++)
		{
			if(i==0) { system[i][k] = 1.0; }
			if(i==1){ system[i][k] = zero[k];}
			else
			{
				system[i][k] = pow(zero[k], i);
			}
		}
	}
	gaussianElimination(system, w, n);
	#pragma omp parallel for shared(m) reduction(+:gqIntegral)		
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			for(k=0;k<n;k++)
			{
				for(l=0;l<n;l++)
				{
					gqIntegral += w[i] * w[j] * w[k] * w[l] * f(zero[l], zero[k], zero[j], zero[i], m);
				}
			}
		}
	}		
	delete[] w, zero, A, a;
	return gqIntegral;
}

int main()
{
	int i, n = 0;	
	double m = 0.0;
	double Integral = 0.0;
	for(m=pow(10.0,-0.5); m<=10.0;m*=pow(10.0,1.0/6.0))
	{
		for(n=2;n<=15;n++)
		{
			Integral = gaussQuad(n, m);
			printf("%d %.15f\n", n, Integral);
		}
	}
	return 0;
}
	