
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <stdlib.h>

// 1D length
#define N 2048

// Maximum number of iterations
#define ITER_MAX 1000000

// How often to check the relative residual
#define RESID_FREQ 1000 

// The residual
#define RESID 1e-6

// Useful globals
int world_size; // number of processes
int my_rank; // my process number

double magnitude(double** x, const int size);
void jacobi(double** x, double** b, double** tmp, const int size);
double getResid(double**x, double** b, const int size);

int main(int argc, char** argv)
{

   int i,j, totiter;
   int done = 0;
   double **x, **xtemp, **b;
   double bmag, resmag;
   int Nrows; 

   // Initialize MPI
   MPI_Init(&argc, &argv);
	   
   // Get the number of processes
   MPI_Comm_size(MPI_COMM_WORLD, &world_size);
   
   // Get the rank
   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
   
   // Figure out my local size. The last rank gets the leftover. 
   Nrows = N/world_size;
   
   if (my_rank == (world_size-1)) { Nrows += (N % world_size) ; }

   double** x = new double*[Nrows];
   double** xtmp = new double*[Nrows];
   double** b = new double*[Nrows];
   for(i=0;i<N<i++)
   {
      x[i] = new double[N];
      xtmp[i] = new double[N];
      b[i] = new double[N];

   }
  
   for (i=0;i<N;i++) 
      { 
         for(j=0;j<N+2;j++)
      {
         x[i][j] = 0.0; 
         xtmp[i][j] = 0.0; 
         b[i][j] = 0.0; 
      }
   }
 
   // The source only lives on a particular rank!
   int source_rank = (N/2)/(N/world_size);

   if (my_rank == source_rank) { b[N/2 - source_rank*Nrows][N/2 - source_rank*Nrows] = 1.0; }

   //Get magnitude of rhs
   bmag = magnitude(b, Nrows);
   printf("bmag: %.8e\n", bmag);

   for (totiter=RESID_FREQ;totiter<ITER_MAX && done==0;totiter+=RESID_FREQ)
   {

      // do RESID_FREQ jacobi iterations
      jacobi(x, b, xtmp, Nrows);

      resmag = getResid(x, b);

      printf("%d res %.8e bmag %.8e rel %.8e\n", totiter, resmag, bmag, resmag/bmag);
      if (resmag/bmag < RESID) { done = 1; }
   }
   
   // Clean up
   MPI_Finalize();

   
   return 0;
}

double magnitude(double** x, const int size)
{
   int i, j;
   double bmag;
   double global_bmag; // used for global reduce! 

   i, j = 0;
   bmag = 0.0;  
   global_bmag = 0.0;
   for (i=0; i<N; i++)
   {
   	for(j=0; j<N;j++)
   	{
     	   bmag = bmag + x[i][j]*x[i][j];
   	}
   }

   // Reduce. 
   MPI_Allreduce(&bmag, &global_bmag, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
   
   return sqrt(global_bmag);
}

void jacobi(double** x, double** b, double** tmp, const int size)
{
   int iter,i, j;

   iter = 0; i = 0;

   for (iter=0;iter<RESID_FREQ;iter++)
   {
      for (i=0;i<N;i++)
      {
      	for(j=0;j<N;j++)
      	{
            tmp[i][j] = (1/4)*(tmp[i+1][j] + tmp[i-1][j] + tmp[i][j+1] + tmp[i][j-1]) + b[i][j];
      	}
      }

      for (i=1;i<=N;i++)
      {
      	for(j=1;j<=N;j++)
      	{
            x[i][j] = tmp[i][j];
      	}
      }
   }
}

double getResid(double** x, double** b)
{
   int i, j;
   double localres,resmag;

   i = 0;
   localres = 0.0;
   resmag = 0.0;

   for (i=1;i<=N;i++)
   {
   	for(j=1;j<=N;j++)
   	{
         localres = (b[i][j] - x[i][j] + 0.5*(x[i+1][j] + x[i-1][j] + x[i][j+1] + x[i][j-1]));
         localres = localres*localres;
         resmag = resmag + localres;
  	}
   }

   resmag = sqrt(resmag);

   return resmag;
}