#include <iostream>
#include <fstream>
#define N 115

using namespace std;


int main()
{
	string line;
	ifstream fin("Complete_TAVG_summary.txt");

	double* fx = new double[N];
	double* gx = new double[N];
	int* x = new int[N];
	double* y = new double[N];

	while(getline(fin,line))
	{
	}

	for(int i=0; i<=115; i++)
	{
		x[i] = 1900 + i;
	}


	for(k=0; k<=N;k++)
	{
		fx[k] = 0;
		for(int j=0;j<N;j++)
		{
			gx[k] = 1.0;
			for(i=0;i<N,i!=j;i++)
			{
				gx[k] *= (x[k]-x[i]) / (x[j] - x[i]);
			}
			fx[k] += y[j] * gx[k];
		}
	}
	
	return 0;
}
