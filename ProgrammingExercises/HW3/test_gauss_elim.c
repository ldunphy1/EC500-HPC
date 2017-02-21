#include <stdio.h>
#include <stdlib.h>

using namespace std;

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
	//delete[] tempA;
	double* temp = new double[dim];
	for(i=0;i<dim;i++)
	{
		temp[i] = 1.0;
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
	//delete[] temp;
}

int main()
{
	int i = 0;
	int dim = 3;
	double** A =  new double*[dim];
	double* b = new double[dim];
	double** AA = new double*[dim];
	double* bb = new double[dim];
	for (i=0;i<dim;i++) { A[i] = new double[dim];}
	A[0][0] = -1.0; A[0][1] = 2.0; A[0][2] = -5.0; b[0] = 17.0;
	A[1][0] = 2.0; A[1][1] = 1.0; A[1][2] = 3.0; b[1] = 0.0;
	A[2][0] = 4.0; A[2][1] = -3.0; A[2][2] = 1.0; b[2] = -10.0;
	gaussianElimination(A,b,dim);
	for (i=0;i<dim;i++) { printf("%f ", b[i]); }
	printf("\n");
	/*delete[] b;
	for (i = 0; i < dim; i++) { delete[] A[i]; }
	delete[] A;*/
	dim = 4;
	for (i=0;i<dim;i++) { AA[i] = new double[dim];}
	AA[0][0] = 1.0; AA[0][1] = 2.0; AA[0][2] = 3.0; AA[0][3] = 4.0; bb[0] = 1.0;
	AA[1][0] = 4.0; AA[1][1] = 8.0; AA[1][2] = 6.0; AA[1][3] = 7.0; bb[1] = 2.0;
	AA[2][0] = 7.0; AA[2][1] = 8.0; AA[2][2] = 10.0; AA[2][3] = 11.0; bb[2] = 3.0;
	AA[3][0] = 1.0; AA[3][1] = 1.0; AA[3][2] = 1.0; AA[3][3] = 5.0; bb[3] = 4.0;
	gaussianElimination(AA,bb,dim);
	for (i=0;i<dim;i++) { printf("%f ", bb[i]); }
	printf("\n");
	/*delete[] bb;
	for (i = 0; i < dim; i++) { delete[] AA[i]; }
	delete[] AA;*/
	
	return 0;
}