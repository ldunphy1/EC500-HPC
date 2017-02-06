#include <iostream>
#include <iomanip> 
#include <cmath>
#include <float.h>

using namespace std;

int getLegendreCoeff(double* A, int n)
{
	double (*a)[n][n] = (double (*)[n][n]) A;
	if((n<0) || (a==NULL))
		return 0;
	(*a)[0][0] = 1;
	(*a)[1][0] = 0;
	(*a)[1][1] = 1;
	for(int j=n;j>=0;j--)
	{
		(*a)[n][j] = ((2*n-1) * (*a)[n-1][j-1] - (n-1) * (*a)[n-2][j])/n;		
	}
	return 1;		
}

int main()
{
	int n;
	cout<<"What order Legendre polynomial? "<<endl;
	cin>>n;
	double a[n+1][n+1];	
	int x = getLegendreCoeff(*a,n);
	cout<<a[n][0]<<" x^"<<n;
	for(int i=n-1;i>=0;i--)
	{
		cout<<" + "<<a[n][i]<<" x^"<<i;
	}
	cout<<endl;
	
}