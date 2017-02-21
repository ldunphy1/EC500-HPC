#include <iostream>
#include <iomanip> 
#include <cmath>
#include <float.h>
#define PI 3.14159265358979323846
#define TOL 1.0E-8

using namespace std;

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

int main()
{
	int n;
	cout<<"What order Legendre polynomial?\n";
	cin>>n;
	if(n>30)
	{
		cout<<"C'mon man, don't ask me for an order greater than 30.\n";
		cin>>n;
	}
	double* zero = new double[n+1];
	double* a = new double[n+1];
	for(int i=0; i<n+1;i++)
	{
		zero[i] = 0.0;
		a[i] = 0.0;
	}
	double* A = new double[(n+1)*(n+1)];
	for(int i=0;i<(n+1)*(n+1);i++)
	{
		A[i] =0.0;
	}
	getLegendreCoeff(A, n);
	for(int i=0;i<n+1;i++)
	{
		a[i] =A[i + n*(n+1)];
	}
	delete[] A;
	getLegendreZero(zero, a, n);
	delete[] a;
	for(int i=0;i<n;i++)
	{
		cout<<zero[i]<<" ";
	}
	cout<<endl;
	delete[] zero;
	
	return 0;
}
