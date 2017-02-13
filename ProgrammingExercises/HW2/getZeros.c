#include <iostream>
#include <iomanip> 
#include <cmath>
#include <float.h>
#include "legendre.c"

using namespace std;

int getLegendreZero(double* zero, double* a, int n)
{
	const double PI = 3.141592653589793;
	int k = (1.0 - (1.0/8.0*pow(n,2.0))+(1.0/(8.0*pow(n,3.0))))*cos(PI*(4*)
	int y = getLegendreCoeff(a, n);
}

int main()
{
	int n;
	cout<<"Wat order Legendre polynomial?\n";
	cin>>n;
	if(n>30)
	{
		cout<<"C'mon man, don't ask me for an order greater than 30.\n";
		cin>>n;
	}
	double* zero = new double(n);
	double* a = new double(n+1);
	for(int i=0; i<n+1;i++)
	{
		zero[i] = 0.0;
		a[i] = 0.0;
	}
	int x = getLegendreZero(zero, a, n);
	
	return 0;
}