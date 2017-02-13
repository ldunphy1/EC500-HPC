#include <stdio.h>
#include <stdlib.h>

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

int main()
{
	int i;
	double** A;
	double* b;
	int dim = 3;
	A = new double*[dim];
	for (i=0;i<dim;i++) { A[i] = new double[dim];}
	b = new double[dim];
	A[0][0] = -1; A[0][1] = 2; A[0][2] = -5; b[0] = 17;
	A[1][0] = 2; A[1][1] = 1; A[1][2] = 3; b[1] = 0;
	A[2][0] = 4; A[2][1] = -3; A[2][2] = 1; b[2] = -10;
	gaussianElimination(A,b,dim);
	for (i=0;i<dim;i++) { printf("%f ", b[i]); }
	printf("\n");
	delete[] b;
	for (i = 0; i < dim; i++) { delete[] A[i]; }
	delete[] A;
	dim = 4;
	A = new double*[dim];
	for (i=0;i<dim;i++) { A[i] = new double[dim];}
	b = new double[dim];
	A[0][0] = 1; A[0][1] = 2; A[0][2] = 3; A[0][3] = 4; b[0] = 1;
	A[1][0] = 4; A[1][1] = 8; A[1][2] = 6; A[1][3] = 7; b[1] = 2;
	A[2][0] = 7; A[2][1] = 8; A[2][2] = 10; A[2][3] = 11; b[2] = 3;
	A[3][0] = 1; A[3][1] = 1; A[3][2] = 1; A[3][3] = 5; b[3] = 4;
	gaussianElimination(A,b,dim);
	for (i=0;i<dim;i++) { printf("%f ", b[i]); }
	printf("\n");
	delete[] b;
	for (i = 0; i < dim; i++) { delete[] A[i]; }
	delete[] A;
	
	return (0);
}