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
	for(int i=0;i<=n;i++)
	{
		(*a)[0][i] = 1;
		(*a)[1][i] = i;
	}
	for(int i=0;i<=n;i++)
	{
		for(int j=0;j<=n;j++)
		{
			(*a)[i][j] = (2*i-1)/i * j * (*a)[i-1][j] - (i-1)/i * (*a)[i-2][j];		
		}
	}
	return 1;		
}

int main()
{
	int n;
	cout<<"What order Legendre polynomial? "<<endl;
	cin>>n;
	double a[n+1][n+1];	
	int x = getLegendreCoeff((*a),n);
	cout<<a[n][0]<<" x^"<<n;
	for(int i=n-1;i>=0;i--)
	{
		cout<<" + "<<a[n][i]<<" x^"<<i;
	}
	cout<<endl;
	
}