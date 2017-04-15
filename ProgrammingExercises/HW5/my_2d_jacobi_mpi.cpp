
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <stdlib.h>

// Maximum number of iterations
#define ITER_MAX 10000000

// How often to check the relative residual
#define RESID_FREQ 10000

// The residual
#define RESID 1e-2

// Useful globals
int world_size; // number of processes
int my_rank;    // my process number
int N;

double magnitude(double **x, const int size);
void jacobi(double **x, double **b, double **tmp, const int size);
double getResid(double **x, double **b, const int size);

int main(int argc, char **argv)
{
      // Initialize MPI
      MPI_Init(&argc, &argv);

      for (N = 16; N <= 512; N += 2)
      {

            int i, j, totiter;
            int done = 0;
            double bmag, resmag;
            int Nrows = 0;

            // Get the number of processes
            MPI_Comm_size(MPI_COMM_WORLD, &world_size);

            // Get the rank
            MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

            // Figure out my local size. The last rank gets the leftover.
            Nrows = N / world_size;

            if (my_rank == (world_size - 1))
            {
                  Nrows += (N % world_size);
            }

            double **x = new double *[Nrows + 1];
            double **xtmp = new double *[Nrows + 1];
            double **b = new double *[Nrows + 1];
            for (i = 0; i < N + 1; i++)
            {
                  x[i] = new double[N + 1];
                  xtmp[i] = new double[N + 1];
                  b[i] = new double[N + 1];
            }

            for (i = 0; i < N + 1; i++)
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
                  b[N / 2 - source_rank * Nrows][N / 2 - source_rank * Nrows] = 1.0;
            }

            //Get magnitude of rhs
            bmag = magnitude(b, Nrows);
            printf("bmag: %.8e\n", bmag);

            for (totiter = RESID_FREQ; totiter < ITER_MAX && done == 0; totiter += RESID_FREQ)
            {
                  printf("can enter loop\n");
                  // do RESID_FREQ jacobi iterations
                  jacobi(x, b, xtmp, Nrows);
                  printf("finished jacobi\n");

                  resmag = getResid(x, b, Nrows);
                  printf("finished getResid\n");

                  printf("%d res %.8e bmag %.8e rel %.8e\n", totiter, resmag, bmag, resmag / bmag);
                  if (resmag / bmag < RESID)
                  {
                        done = 1;
                  }
            }
            for (i = 0; i < N + 1; i++)
            {
                  delete[] x[i], xtmp[i], b[i];
            }
            delete[] x, xtmp, b;
      }

      // Clean up
      MPI_Finalize();

      return 0;
}

double magnitude(double **x, const int size)
{
      int i;
      int j;
      double bmag;
      double global_bmag; // used for global reduce!
      const int lower_limit = (my_rank == 0) ? 1 : 0;
      const int upper_limit = (my_rank == world_size - 1) ? size - 1 : size;

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

void jacobi(double **x, double **b, double **tmp, const int size)
{
      printf("entered jacobi\n");
      int iter, i, j;

      // Prepare for async send/recv
      MPI_Request request[4];
      int requests;
      MPI_Status status[4];

      const int lower_limit = (my_rank == 0) ? 1 : 0;
      const int upper_limit = (my_rank == world_size - 1) ? size - 1 : size;

      // grab the left and right buffer.
      double left_buffer = 0.0;
      double right_buffer = 0.0;

      iter = 0;
      i = 0;

      for (iter = 0; iter < RESID_FREQ; iter++)
      {
            requests = 0;
            // Fill the left buffer. Send to the right, listen from the left.
            MPI_Isend(&x[size-1][0], N + 1, MPI_DOUBLE, (my_rank + 1) % world_size, 1, MPI_COMM_WORLD, request + requests++);
            MPI_Irecv(&left_buffer, 1, MPI_DOUBLE, (my_rank + world_size - 1) % world_size, 1, MPI_COMM_WORLD, request + requests++);

            printf("send & receive 1\n");

            // Fill the right buffer. Send to the left, listen from the right.
            MPI_Isend(&x[0][0], N + 1, MPI_DOUBLE, (my_rank + world_size - 1) % world_size, 0, MPI_COMM_WORLD, request + requests++);
            MPI_Irecv(&right_buffer, 1, MPI_DOUBLE, (my_rank + 1) % world_size, 0, MPI_COMM_WORLD, request + requests++);

            printf("send & receive 2\n");

            for (i = 1; i < size; i++)
            {
                  for (j = 1; j < N; j++)
                  {
                        tmp[i][j] = (1 / 4) * (tmp[i + 1][j] + tmp[i - 1][j] + tmp[i][j + 1] + tmp[i][j - 1]) + b[i][j];
                  }
            }

            printf("before wait all\n");

            // Wait for async.
            MPI_Waitall(requests, request, status);

            printf("after wait all\n");

            // Impose zero bc.
            if (my_rank != 0)
            {
                  for (i = 0; i < N + 1; i++)
                  {
                        tmp[0][i] = 0.5 * (x[0][i] + left_buffer) + b[0][i];
                  }
            }

            // Impose zero bc.
            if (my_rank != world_size - 1)
            {
                  for (i = 0; i < N + 1; i++)
                  {
                        tmp[size - 1][i] = 0.5 * (right_buffer + x[size - 1][i]) + b[size - 1][i];
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

      printf("before barrier\n");

      MPI_Barrier(MPI_COMM_WORLD);

      printf("after barrier\n");
}

double getResid(double **x, double **b, const int size)
{
      const int lower_limit = (my_rank == 0) ? 1 : 0;
      const int upper_limit = (my_rank == world_size - 1) ? size - 1 : size;

      int i, j;
      double localres, resmag;
      double global_resmag;

      // Prepare for async send/recv
      MPI_Request request[4];
      int requests;
      MPI_Status status[4];

      // grab the left and right buffer.
      double left_buffer = 0.0;
      double right_buffer = 0.0;

      requests = 0;

      // Fill the left buffer. Send to the right, listen from the left.
      MPI_Isend(&x[size-1][0], N + 1, MPI_DOUBLE, (my_rank + 1) % world_size, 1, MPI_COMM_WORLD, request + requests++);
      MPI_Irecv(&left_buffer, 1, MPI_DOUBLE, (my_rank + world_size - 1) % world_size, 1, MPI_COMM_WORLD, request + requests++);

      // Fill the right buffer. Send to the left, listen from the right.
      MPI_Isend(&x[0][0], N + 1, MPI_DOUBLE, (my_rank + world_size - 1) % world_size, 0, MPI_COMM_WORLD, request + requests++);
      MPI_Irecv(&right_buffer, 1, MPI_DOUBLE, (my_rank + 1) % world_size, 0, MPI_COMM_WORLD, request + requests++);

      i = 0;
      localres = 0.0;
      resmag = 0.0;

      for (i = lower_limit; i < upper_limit; i++)
      {
            for (j = 1; j <= N; j++)
            {
                  localres = (b[i][j] - x[i][j] + 0.5 * (x[i + 1][j] + x[i - 1][j] + x[i][j + 1] + x[i][j - 1]));
                  localres = localres * localres;
                  resmag = resmag + localres;
            }
      }

      // Wait for async.
      MPI_Waitall(requests, request, status);

      // Impose zero bc.
      if (my_rank != 0)
      {
            for (i = 0; i < N + 1; i++)
            {
                  localres = (b[0][i] - x[0][i] + .5 * (x[0][i] + left_buffer));
                  localres = localres * localres;
                  resmag = resmag + localres;
            }
      }

      // Impose zero bc.
      if (my_rank != world_size - 1)
      {
            for (i = 0; i < N + 1; i++)
            {
                  localres = (b[size - 1][i] - x[size - 1][i] + 0.5 * (right_buffer + x[size - 1][i]));
                  localres = localres * localres;
                  resmag = resmag + localres;
            }
      }

      // Reduce.
      MPI_Allreduce(&resmag, &global_resmag, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      return sqrt(global_resmag);
}