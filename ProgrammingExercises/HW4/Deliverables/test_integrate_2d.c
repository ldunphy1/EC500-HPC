#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#define PI 3.14159265358979323846
#define TOL 1.0E-8

using namespace std;
double g1(double x, double y)
{
	return (pow(x, 8.0) + pow(y, 8.0) + pow((y - 1.0), 3.0) * pow((x - 3.0), 5.0));
}

double g2(double x, double y)
{
	return sqrt(pow(x, 2.0) - pow(y, 2.0) + 2.0);
}

double g3(double x, double y)
{
	return (exp(-pow(x, 2.0) - (pow(y, 2.0)/8.0)) * cos(PI * x) * sin(x * PI/8.0));
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

void trapRule(int n, double* trapIntegrals)
{
	int i = 0;
	int j = 0;
	double h = 2.0/n;
	double trapInt1A = 0.0;
	double trapInt1B = 0.0;
	double trapInt1C = 0.0;
	double trapInt2A = 0.0;
	double trapInt2B = 0.0;
	double trapInt2C = 0.0;
	double trapInt3A = 0.0;
	double trapInt3B = 0.0;
	double trapInt3C = 0.0;
	double* x = new double[n+1];
	double* y = new double[n+1];
	for (i = 0; i < 3; i++)
	{
		trapIntegrals[i] = 0.0;
	}		
	x[0] = -1.0;
	x[n] = 1.0; 
	y[0] = -1.0;
	y[n] = 1.0;
	for(i=1;i<n;i++)
	{
		x[i] = -1.0 + (2.0 * i/n);
		y[i] = -1.0 + (2.0 * i/n);
	}
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			trapInt1A += g1(x[j], y[i]);
			trapInt2A += g2(x[j], y[i]);
			trapInt3A += g3(x[j], y[i]);
			trapInt1B += g1(x[0], y[i]) + g1(x[n], y[i]);
			trapInt2B += g2(x[0], y[i]) + g2(x[n], y[i]);
			trapInt3B += g3(x[0], y[i]) + g3(x[n], y[i]);
			trapInt1C += g1(x[j], y[0]) + g1(x[j], y[n]);
			trapInt2C += g2(x[j], y[0]) + g2(x[j], y[n]);
			trapInt3C += g3(x[j], y[0]) + g3(x[j], y[n]);
		}
	}
	trapIntegrals[0] = (pow(h, 2.0)/4) * (g1(x[0], y[0]) + g1(x[0], y[n]) + g1(x[n], y[n]) + 4*trapInt1A + 2*trapInt1B + 2*trapInt1C);
	trapIntegrals[1] = (pow(h, 2.0)/4) * (g2(x[0], y[0]) + g2(x[0], y[n]) + g2(x[n], y[n]) + 4*trapInt2A + 2*trapInt2B + 2*trapInt2C);
	trapIntegrals[2] = (pow(h, 2.0)/4) * (g3(x[0], y[0]) + g3(x[0], y[n]) + g3(x[n], y[n]) + 4*trapInt3A + 2*trapInt3B + 2*trapInt3C);
	delete[] x, y;
}
void gaussQuad(int n, double* gqIntegrals)
{
	int i, j, k = 0;
	for (i = 0; i < 3; i++)
	{
		gqIntegrals[i] = 0.0;
	}
	double* w = new double[n];
	double** system = new double*[n];
	double* zero = new double[n+1];
	double* a = new double[n+1];
	double* A = new double[(n+1)*(n+1)];
	for(int i=0; i<n+1;i++)
	{
		zero[i] = 0.0;
		a[i] = 0.0;
	}
	for(int i=0;i<(n+1)*(n+1);i++)
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
	for(j=0;j<n;j++)
	{
		for(i=0; i<n; i++)
		{
			gqIntegrals[0] += w[i] * w[j] * g1(zero[i], zero[j]);
			gqIntegrals[1] += w[i] * w[j] * g2(zero[i], zero[j]);
			gqIntegrals[2] += w[i] * w[j] * g3(zero[i], zero[j]);
		}
	}
	delete[] w, zero, A, a;
}

int main()
{
	int i, n = 0;	
	double* gqIntegrals = new double[3];
	double* trapIntegrals = new double[3];
	// double* exact = new double[3];
	// double* error = new double[3];
	// exact[0] = 2688.88888888889
	// exact[1] = 5.62421;
	// exact[2] = 0;
	for (i = 0; i < 3; i++)
	{
		gqIntegrals[i] = 0.0;
		trapIntegrals[i] = 0.0;
	}
	printf("approximating using the Trapezoidal Rule\n");
	printf("%s %s %25s %15s\n", "N", "g1", "g2", "g3");
	for(n=2;n<=15;n++)
	{
		trapRule(n, trapIntegrals);
		printf("%d %.15f %.15f %.15f\n", n, trapIntegrals[0], trapIntegrals[1], trapIntegrals[2]);
		//gaussQuad(n, gqIntegrals);
		//printf("%d %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f\n", n, trapIntegrals[0], gqIntegrals[0], exact[0], trapIntegrals[1], gqIntegrals[1], exact[1], trapIntegrals[2], gqIntegrals[2], exact[2]);
	}
	printf("\napproximating using Gaussian Quadrature\n");
	printf("%s %s %25s %15s\n", "N", "g1", "g2", "g3");
	for(n=2;n<=15;n++)
	{
		gaussQuad(n, gqIntegrals);
		printf("%d %.15f %.15f %.15f\n", n, gqIntegrals[0], gqIntegrals[1], gqIntegrals[2]);
	}	
	delete[] gqIntegrals, trapIntegrals;

	return 0;
}
	