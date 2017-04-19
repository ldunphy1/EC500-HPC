
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <stdlib.h>

// Maximum number of iterations
#define ITER_MAX 10000000

// How often to check the relative residual
#define RESID_FREQ 1000

// The residual
#define RESID 1e-6

// Useful globals
int world_size; // number of processes
int my_rank;    // my process number

double magnitude(double **x, const int Nrows, const int N);
void jacobi(double **x, double **b, double **tmp, const int Nrows, const int N);
double getResid(double **x, double **b, const int Nrows, const int N);

int main(int argc, char **argv)
{
      // Initialize MPI
      MPI_Init(&argc, &argv);

      // Get the number of processes
      MPI_Comm_size(MPI_COMM_WORLD, &world_size);

      // Get the rank
      MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

      if (my_rank == 0)
      {
            printf("I am rank %d of %d, and I speak for the trees!\n", my_rank, world_size);
            fflush(stdout);
      }

      for (int N = 16; N <= 512; N *= 2)
      {
            int i, j, totiter;
            int done = 0;
            double bmag, resmag;
            int Nrows = 0;
            
            if (my_rank == 0)
            {
                  printf("We're looking at N = %d.\n", N);
                  fflush(stdout);
            }

            // Figure out my local size. The last rank gets the leftover.
            Nrows = N / world_size;

            if (my_rank == (world_size - 1))
            {
                  Nrows += (N % world_size);
            }

            double **x = new double *[Nrows];
            double **xtmp = new double *[Nrows];
            double **b = new double *[Nrows];
            for (i = 0; i < Nrows; i++)
            {
                  x[i] = new double[N + 1];
                  xtmp[i] = new double[N + 1];
                  b[i] = new double[N + 1];
            }

            for (i = 0; i < Nrows; i++)
            {
                  for (j = 0; j < N + 1; j++)
                  {
                        x[i][j] = 0.0;
                        xtmp[i][j] = 0.0;
                        b[i][j] = 0.0;
                  }
            }

            // The source only lives on a particular rank!
            int source_rank = (N / 2) / (N / world_size);

            if (my_rank == source_rank)
            {
                  b[N / 2 - source_rank * Nrows][N / 2] = 1.0;
            }

            //Get magnitude of rhs
            bmag = magnitude(b, Nrows, N);
            if (my_rank == 0)
            {
                  printf("bmag: %.8e\n", bmag);
                  fflush(stdout);
            }

            for (totiter = RESID_FREQ; totiter < ITER_MAX && done == 0; totiter += RESID_FREQ)
            {
                  //if (my_rank == 0) printf("can enter loop\n");
                  // do RESID_FREQ jacobi iterations
                  jacobi(x, b, xtmp, Nrows, N);
                  //if (my_rank == 0) printf("finished jacobi\n");

                  resmag = getResid(x, b, Nrows, N);
                  //if (my_rank == 0) printf("finished getResid\n");

                  if (my_rank == 0)
                  {
                        printf("%d res %.8e bmag %.8e rel %.8e\n", totiter, resmag, bmag, resmag / bmag);
                        fflush(stdout);
                  }
                  if (resmag / bmag < RESID)
                  {
                        done = 1;
                  }
            }
            for (i = 0; i < Nrows; i++)
            {
                  delete[] x[i], xtmp[i], b[i];
            }
            delete[] x, xtmp, b;
      }

      // Clean up
      MPI_Finalize();

      return 0;
}

double magnitude(double **x, const int Nrows, const int N)
{
      int i;
      int j;
      double bmag;
      double global_bmag; // used for global reduce!
      const int lower_limit = (my_rank == 0) ? 1 : 0;
      const int upper_limit = (my_rank == world_size - 1) ? Nrows - 1 : Nrows;

      i, j = 0;
      bmag = 0.0;
      global_bmag = 0.0;
      for (i = lower_limit; i < upper_limit; i++)
      {
            for (j = 1; j < N; j++)
            {
                  bmag = bmag + x[i][j] * x[i][j];
            }
      }

      // Reduce.
      MPI_Allreduce(&bmag, &global_bmag, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      return sqrt(global_bmag);
}

void jacobi(double **x, double **b, double **tmp, const int Nrows, const int N)
{
      int iter, i, j;

      // Prepare for async send/recv
      MPI_Request request[4];
      int requests;
      MPI_Status status[4];

      const int lower_limit = (my_rank == 0) ? 1 : 0;
      const int upper_limit = (my_rank == world_size - 1) ? Nrows - 1 : Nrows;

      // grab the top and bottom buffer.
      double* top_buffer = new double[N + 1];
      double* bottom_buffer = new double[N + 1];

      iter = 0;
      i = 0;

      for (iter = 0; iter < RESID_FREQ; iter++)
      {
            requests = 0;
            // Fill the top buffer. Send to the bottom, listen from the top
            MPI_Isend(x[Nrows-1], N + 1, MPI_DOUBLE, (my_rank + 1) % world_size, 1, MPI_COMM_WORLD, request + requests++);
            MPI_Irecv(top_buffer, N + 1, MPI_DOUBLE, (my_rank + world_size - 1) % world_size, 1, MPI_COMM_WORLD, request + requests++);


            // Fill the bottom buffer. Send to the top, listen from the bottom. 
            MPI_Isend(x[0], N + 1, MPI_DOUBLE, (my_rank + world_size - 1) % world_size, 0, MPI_COMM_WORLD, request + requests++);
            MPI_Irecv(bottom_buffer, N + 1, MPI_DOUBLE, (my_rank + 1) % world_size, 0, MPI_COMM_WORLD, request + requests++);

            for (i = 1; i < Nrows-1; i++)
            {
                  for (j = 1; j < N; j++)
                  {
                        tmp[i][j] = 0.25 * (x[i + 1][j] + x[i - 1][j] + x[i][j + 1] + x[i][j - 1]) + b[i][j];
                  }
            }

            // Wait for async.
            MPI_Waitall(requests, request, status);

            // Impose zero bc.
            if (my_rank != 0)
            {
                  for (j = 1; j < N; j++)
                  {
                        tmp[0][j] = 0.25 * (x[1][j] + top_buffer[j] + x[0][j + 1] + x[0][j - 1]) + b[0][j];
                  }
            }

            // Impose zero bc.
            if (my_rank != world_size - 1)
            {
                  for (j = 1; j < N; j++)
                  {
                        tmp[Nrows-1][j] = 0.25 * (x[Nrows-2][j] + bottom_buffer[j] + x[Nrows-1][j + 1] + x[Nrows-1][j - 1]) + b[Nrows-1][j];
                  }
            }

            for (i = lower_limit; i < upper_limit; i++)
            {
                  for (j = 1; j < N; j++)
                  {
                        x[i][j] = tmp[i][j];
                  }
            }
      }

      MPI_Barrier(MPI_COMM_WORLD);

      delete[] top_buffer;
      delete[] bottom_buffer;
}

double getResid(double **x, double **b, const int Nrows, const int N)
{
      int i, j;
      double localres, resmag;
      double global_resmag;

      // Prepare for async send/recv
      MPI_Request request[4];
      int requests;
      MPI_Status status[4];

      const int lower_limit = (my_rank == 0) ? 1 : 0;
      const int upper_limit = (my_rank == world_size - 1) ? Nrows - 1 : Nrows;

      // grab the top and bottom buffer.
      double* top_buffer = new double[N + 1];
      double* bottom_buffer = new double[N + 1];

      requests = 0;
      // Fill the top buffer. Send to the bottom, listen from the top
      MPI_Isend(x[Nrows-1], N + 1, MPI_DOUBLE, (my_rank + 1) % world_size, 1, MPI_COMM_WORLD, request + requests++);
      MPI_Irecv(top_buffer, N + 1, MPI_DOUBLE, (my_rank + world_size - 1) % world_size, 1, MPI_COMM_WORLD, request + requests++);

      // Fill the bottom buffer. Send to the top, listen from the bottom. 
      MPI_Isend(x[0], N + 1, MPI_DOUBLE, (my_rank + world_size - 1) % world_size, 0, MPI_COMM_WORLD, request + requests++);
      MPI_Irecv(bottom_buffer, N + 1, MPI_DOUBLE, (my_rank + 1) % world_size, 0, MPI_COMM_WORLD, request + requests++);

      i = 0;
      localres = 0.0;
      resmag = 0.0;

      for (i = 1; i < Nrows-1; i++)
      {
            for (j = 1; j < N; j++)
            {
                  localres = x[i][j] - (0.25 * (x[i + 1][j] + x[i - 1][j] + x[i][j + 1] + x[i][j - 1]) + b[i][j]);
                  localres = localres * localres;
                  resmag = resmag + localres;
            }
      }


      // Wait for async.
      MPI_Waitall(requests, request, status);

      // Impose zero bc.
      if (my_rank != 0)
      {
            for (j = 1; j < N; j++)
            {
                  localres = x[0][j] - (0.25 * (x[1][j] + top_buffer[j] + x[0][j + 1] + x[0][j - 1]) + b[0][j]);
                  localres = localres * localres;
                  resmag = resmag + localres;
            }
      }

      // Impose zero bc.
      if (my_rank != world_size - 1)
      {
            for (j = 1; j < N; j++)
            {
                  localres = x[Nrows-1][j] - (0.25 * (x[Nrows-2][j] + bottom_buffer[j] + x[Nrows-1][j + 1] + x[Nrows-1][j - 1]) + b[Nrows-1][j]);
            }
      }

      // Reduce.
      MPI_Allreduce(&resmag, &global_resmag, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      delete[] top_buffer;
      delete[] bottom_buffer; 

      return sqrt(global_resmag);
}