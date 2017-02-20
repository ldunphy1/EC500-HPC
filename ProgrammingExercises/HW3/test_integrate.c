#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#define PI 3.14159265358979323846
#define TOL 1.0E-8

using namespace std;

int gaussianElimination(double** A, double* b, int dim)
{ 
	int max, i, j, k;
	double* tempA;
	double tempB, factor, sum;
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
	//delete[] tempA;
	double* temp = new double[dim];
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
	delete[] temp;
    return(0);
}

int getLegendreCoeff(double* a, int n)
{		
	a[0+0*(n+1)] = 1.0; //set P0
	a[1+1*(n+1)] = 1.0; //set P1
    for(int pn=2; pn<n+1;pn++)
    { 
    	a[0 +pn*(n+1)] = -((pn-1)*a[0 + (pn-2)*(n+1)])/(double) pn;
        for(int i =1;i < pn+1;i++)
    		a[i+pn*(n+1)] = ((2*pn-1)*a[i-1 +(pn-1)*(n+1)] - (pn-1)*a[i + (pn-2)*(n+1)])/(double) pn;
    }
	return 1;		
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
		sum += a[i] * i * pow(x,i-1);
	}
	return sum;
}

int getLegendreZero(double* zero, double* a, int n)
{
	int k, count=0;
	double xprev, xnext, f, fp;
	bool initial=true;
	for (k=1;k<=n;k++)
	{
		initial = true;
		xprev = (1.0 - (1.0/(8.0 * pow(n,2.0))) + (1.0/(8.0 * pow(n,3)))) * cos(PI * (4.0*k - 1.0)/(4.0*n + 2.0));
		xnext = 2.0;
		while ( (abs((xnext-xprev)/xprev) >= TOL) && (xprev!=0) )
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
	return 0;
}

double f1(double x)
{
	return pow(x, 8.0);
}

double f2(double x)
{
	return cos(PI * x/2.0);
}

double f3(double x)
{
	return (1.0/(pow(x,2.0)+1.0));
}
double* trapRule(int n)
{
	int i;
	double h;
	double* trapIntegrals = new double[3];
	h = 2.0/n;
	double* x = new double[n+1];
	x[0] = -1.0;
	x[n] = 1.0; 
	for(i=1;i<n;i++)
	{
		x[i] = -1.0 + (2.0*i/n);
	}
	for(i=0;i<n;i++)
	{
		trapIntegrals[0] += h * (f1(x[i]) + f1(x[i+1]))/2.0;
		trapIntegrals[1] += h * (f2(x[i]) + f2(x[i+1]))/2.0;
		trapIntegrals[2] += h * (f3(x[i]) + f3(x[i+1]))/2.0;
	}

	return trapIntegrals;	
}
double* gaussQuad(int n)
{
	int i, k;
	double* gqIntegrals = new double[3];
	double* w = new double[n+1];
	double** system = new double*[n];
	double* zero = new double[n+1];
	double* a = new double[n+1];
	double* A = new double[(n+1)*(n+1)];
	for (i=0;i<n;i++) { system[i] = new double[n]; }
	getLegendreCoeff(A, n);
	for (i=0;i<n+1;i++) { a[i] =A[i + n*(n+1)]; }
	delete[] A;
	getLegendreZero(zero, a, n);
	delete[] a;
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
			if(i==1){ system[i][k] = zero[k] ;}
			else
			{
				system[i][k] = pow(zero[k], i);
			}
		}
	}
	gaussianElimination(system, w, n);
	//for (i = 0; i < n; i++) { delete[] system[i]; }
	//delete[] system;
	for(i=1;i<=n;i++)
	{
		gqIntegrals[0] += w[i] * f1(zero[i]);
		gqIntegrals[1] += w[i] * f2(zero[i]);
		gqIntegrals[2] += w[i] * f3(zero[i]);
	}
	delete[] w, zero;

	return gqIntegrals;
}

int main()
{
	int n;	
	double* gqIntegrals = new double[3];
	double* trapIntegrals = new double[3];
	printf("approximating using the Trapezoidal Rule\n");
	printf("%s %s %25s %15s\n", "N", "x^8", "cos(pi*x/2)", "1/(x^2+1)");
	for(n=0;n<=10;n++)
	{
		trapIntegrals = trapRule(n);
		printf("%d %.15f %.15f %.15f\n", n, trapIntegrals[0], trapIntegrals[1], trapIntegrals[2]);
	}
	delete[] trapIntegrals;
	printf("\napproximating using Gaussian Quadrature\n");
	printf("%s %s %25s %15s\n", "N", "x^8", "cos(pi*x/2)", "1/(x^2+1)");
	for(n=0;n<=10;n++)
	{
		gqIntegrals = gaussQuad(n);
		printf("%d %.15f %.15f %.15f\n", n, gqIntegrals[0], gqIntegrals[1], gqIntegrals[2]);

	}	

	return 0;
}
	