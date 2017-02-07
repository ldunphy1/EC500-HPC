#include <iostream>
#include <iomanip> 
#include <cmath>
#include <float.h>

using namespace std;

int getLegendreCoeff(double* a, int n)
{
	for(int i=0;i<=n;i++)
	{
		for(int j=0;j<=n;j++)
		{
			a[i+j*(n+1)] = 0.0;
		}
	}	
	a[0+0*(n+1)] = 1.0;
	if(0<n)
		a[1+1*(n+1)] = 1.0;
	for(int i=2;i<=n;i++)
	{
		for(int j=0;j<=i-2;j++)
		{
			a[i+j*(n+1)] = (double)(-i+1)*a[i-2+j*(n+1)]/(double)i;			
		}
		for(int j=1;j<=i;j++)
		{
			a[i+j*(n+1)] = a[i+j*(n+1)]+(double)(i+i-1)*a[i-1+(j-1)*(n+1)]/(double)i;
		}
	}
	return 1;		
}

int main()
{
	int n;
	cout<<"What order Legendre polynomial? "<<endl;
	cin>>n;
	double a[(n+1)*(n+1)];	
	int x = getLegendreCoeff(a,n);
	cout<<a[n+0*(n+1)]<<" x^"<<n;
	for(int i=n-1;i>=0;i--)
	{
		cout<<" + "<<a[n+i*(n+1)]<<" x^"<<i;
	}
	cout<<endl;
	
}