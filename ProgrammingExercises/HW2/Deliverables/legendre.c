#include <iostream>
#include <iomanip> 
#include <cmath>
#include <float.h>

using namespace std;

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

int main()
{
	int n;
	cout<<"What order Legendre polynomial? "<<endl;
	cin>>n;
	double* a = new double((n+1)*(n+1));	
	for(int i=0;i<(n+1)*(n+1);i++)
	{
		a[i] =0.0;
	}
	int x = getLegendreCoeff(a,n);
	cout<<a[n+n*(n+1)]<<" x^"<<n;
	for(int i=n-1;i>=0;i--)
	{
		cout<<" + "<<a[i +  n*(n+1)]<<" x^"<<i;
	}
	cout<<endl;
	
}
