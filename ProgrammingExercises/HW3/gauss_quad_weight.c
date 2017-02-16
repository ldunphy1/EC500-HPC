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
	delete[] tempA;
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

int main()
{
	int n, i, k;
	cout<<"What value of N? \n";
	cin>>n;
	double* x = new double[n];
	double* w = new double[n];
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
			x[k] = 0.0;
			w[k] = x[k];
		}
		else
		{
			x[k] = 2.0/(k+1.0);
			w[k] = x[k];
		}
	}
	for(i=0;i<n;i++)
	{
		for(k=0;k<n;k++)
		{
			if(i==0) { system[i][k] = 1.0; }
			if(i==1){ system[i][k] = x[k] ;}
			else
			{
				system[i][k] = pow(x[k], i);
			}
		}
	}
	gaussianElimination(system, w, n);
	for (i=0;i<n;i++) 
	{
		cout<<"The zeros and weights are given by: \n";
		cout<<"x_i\tw_i\n";
		for(i=0;i<n;i++)
		{
			cout<<x[i]<<" "<<w[i]<<endl; 
		}
	}
	delete[] x, w;
	for (i = 0; i < n; i++) { delete[] system[i]; }
	delete[] system;


	return 0;
}