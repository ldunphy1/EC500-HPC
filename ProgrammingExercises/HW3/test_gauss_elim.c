#include <stdio.h>

using namespace std;

int gaussianElimination(double** a, double* b, int dim)
{ 
    for(j=1; j<=n; j++) /* loop for the generation of upper triangular matrix*/
    {
        for(i=1; i<=n; i++)
        {
            if(i>j)
            {
                c=A[i][j]/A[j][j];
                for(k=1; k<=n+1; k++)
                {
                    A[i][k]=A[i][k]-c*A[j][k];
                }
            }
        }
    }
    x[n]=A[n][n+1]/A[n][n];
    /* this loop is for backward substitution*/
    for(i=n-1; i>=1; i--)
    {
        sum=0;
        for(j=i+1; j<=n; j++)
        {
            sum=sum+A[i][j]*x[j];
        }
        x[i]=(A[i][n+1]-sum)/A[i][i];
    }
    printf("\nThe solution is: \n");
    for(i=1; i<=n; i++)
    {
        printf("\nx%d=%f\t",i,x[i]); /* x1, x2, x3 are the required solutions*/
    }
    return(0);
}

int main()
{
	/*
	   int i,j,k,n;
    float A[20][20],c,x[10],sum=0.0;
    printf("\nEnter the order of matrix: ");
    scanf("%d",&n);
    printf("\nEnter the elements of augmented matrix row-wise:\n\n");
    for(i=1; i<=n; i++)
    {
        for(j=1; j<=(n+1); j++)
        {
            printf("A[%d][%d] : ", i,j);
            scanf("%f",&A[i][j]);
        }
    }
	*/
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
	/*gaussian elimination(A,b,dim);
	for (i=0;i<dim;i++) { printf("%f ", b[i]); }
	delete[] b;
	for (i = 0; i < dim; i++) { delete[] A[i]; }
	delete[] A;
	
	return 0;
}*/