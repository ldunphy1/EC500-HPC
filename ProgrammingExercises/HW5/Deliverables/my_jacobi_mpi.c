// Based on the tutorial from mpitutorial.com/tutorials/mpi-hello-world/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define GIG 1000000000


// Maximum number of iterations
#define ITER_MAX 10000000

// How often to check the relative residual
#define RESID_FREQ 10000

// The residual
#define RESID 1e-2

// Useful globals
int world_size; // number of processes
int my_rank; // my process number

// Timing structures.
timespec time1, time2, timediff; 

// Timing storage.
double timing=0.0;

timespec diff(timespec start, timespec end)
{
  timespec temp;
  if ((end.tv_nsec-start.tv_nsec)<0)
  {
    temp.tv_sec = end.tv_sec-start.tv_sec-1;
    temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
  } else {
    temp.tv_sec = end.tv_sec-start.tv_sec;
    temp.tv_nsec = end.tv_nsec-start.tv_nsec;
  }
  return temp;
}

double magnitude(double* x, const int size);
void jacobi(double* x, double* b, double* tmp, const int size);
double getResid(double* x, double* b, const int size);

int main(int argc, char** argv)
{
	// Initialize MPI
   MPI_Init(&argc, &argv);

	for(int N=16; N<=16384; N*=4)
	{
	   int i,totiter;
	   int done = 0;
	   double *x, *xtmp, *b; 
	   double bmag, resmag;
	   int local_size; 
	   
	   // Get the number of processes
	   MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	   
	   // Get the rank
	   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	   
	   // Figure out my local size. The last rank gets the leftover. 
	   local_size = N/world_size;
	   
	   if (my_rank == (world_size-1)) { local_size += (N % world_size) ; }

	   x = (double*)malloc(local_size*sizeof(double));
	   xtmp = (double*)malloc(local_size*sizeof(double));
	   b = (double*)malloc(local_size*sizeof(double));

	   for (i=0;i<local_size;i++) { x[i] = 0.0; xtmp[i] = 0.0; b[i] = 0.0; }
	   
	   // b[N/2] = 1.0;
	   // The source only lives on a particular rank!
	   int source_rank = (N/2)/(N/world_size);
		//printf("source rank set %d\n", N);
	   //printf("The source at %d goes on rank %d.\n", N/2, source_rank);
	   if (my_rank == source_rank) { b[N/2 - source_rank*(N/world_size)] = 1.0; }
	   //printf("checked my rank=source rank %d\n", N);
	   
	   //Get magnitude of rhs
	   bmag = magnitude(b, local_size);
	   //printf("bmag: %.8e\n", bmag);

	   //printf("I am rank %d of %d and I see the magnitude is %.8e.\n", my_rank, world_size, bmag);



	   for (totiter=RESID_FREQ;totiter<ITER_MAX && done==0;totiter+=RESID_FREQ)
	   {
          // Start timer.
         clock_gettime(CLOCK_REALTIME, &time1);

	      // do RESID_FREQ jacobi iterations
	      jacobi(x, b, xtmp, local_size);        

	      resmag = getResid(x, b, local_size);
	      
	      if (my_rank == 0) {
	         //printf("N %d iterations %d res %.8e bmag %.8e rel %.8e\n", N, totiter, resmag, bmag, resmag/bmag);
	      }
	      
	      if (resmag/bmag < RESID) 
         { 
             // End timer.
            clock_gettime(CLOCK_REALTIME, &time2);

            // Get difference.
            timediff = diff(time1, time2);

            // Save time in seconds to vector.
            timing = ((double)GIG * (double)timediff.tv_sec + (double)timediff.tv_nsec)/((double)GIG);

            done = 1; 
         }
	   }
      printf("%d %.8e\n", N, timing);
	   free(x); free(xtmp); free(b);
	}

	   // Clean up
	   MPI_Finalize();
   
   return 0;
}

double magnitude(double* x, const int size)
{
   int i;
   double bmag;
   double global_bmag; // used for global reduce!
   const int lower_limit = (my_rank == 0) ? 1 : 0;
   const int upper_limit = (my_rank == world_size-1) ? size-1 : size;

   i = 0;
   bmag = 0.0;  
   global_bmag = 0.0;
   for (i = lower_limit; i<upper_limit; i++)
   {
     bmag = bmag + x[i]*x[i];
   }
   
   //printf("I am rank %d of %d and I have a local magnitude %.8e.\n", my_rank, world_size, bmag);

   // Reduce. 
   MPI_Allreduce(&bmag, &global_bmag, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
   
   return sqrt(global_bmag);
}

void jacobi(double* x, double* b, double* tmp, const int size)
{
   int iter,i;
   
   // Prepare for async send/recv
   MPI_Request request[4];
   int requests;
   MPI_Status status[4];
   
   const int lower_limit = (my_rank == 0) ? 1 : 0;
   const int upper_limit = (my_rank == world_size-1) ? size-1 : size;
   
   // grab the left and right buffer.
   double left_buffer = 0.0;
   double right_buffer = 0.0;
   

   iter = 0; i = 0;

   {
      for (iter=0;iter<RESID_FREQ;iter++)
      {

         requests=0;
      
         // Fill the left buffer. Send to the right, listen from the left.
         MPI_Isend(&x[size-1],   1, MPI_DOUBLE, (my_rank+1)%world_size, 1, MPI_COMM_WORLD, request + requests++);
         MPI_Irecv(&left_buffer, 1, MPI_DOUBLE, (my_rank+world_size-1)%world_size, 1, MPI_COMM_WORLD, request + requests++);

         //printf("I am rank %d of %d and I received %.8e from the left.\n", my_rank, world_size, left_buffer);
         
         // Fill the right buffer. Send to the left, listen from the right.
         MPI_Isend(&x[0],   1, MPI_DOUBLE, (my_rank+world_size-1)%world_size, 0, MPI_COMM_WORLD, request + requests++);
         MPI_Irecv(&right_buffer, 1, MPI_DOUBLE, (my_rank+1)%world_size, 0, MPI_COMM_WORLD, request + requests++);
         //printf("I am rank %d of %d and I received %.8e from the right.\n", my_rank, world_size, right_buffer);
         
        
         // Loop over the rest.
         for (i=1;i<size-1;i++)
         {
            tmp[i] = 0.5*(x[i+1]+x[i-1]) + b[i];
         }

         // Wait for async.
         MPI_Waitall ( requests, request, status );
         
         
         // Impose zero bc.
         if (my_rank != 0)
         {
            tmp[0] = 0.5*(x[1]+left_buffer) + b[0];
         }
         
         // Impose zero bc.
         if (my_rank != world_size-1)
         {
            tmp[size-1] = 0.5*(right_buffer+x[size-2]) + b[size-1];
         }
          

         for (i=lower_limit;i<upper_limit;i++)
         {
            x[i] = tmp[i];
         }
         
      }
   }
   
   MPI_Barrier(MPI_COMM_WORLD);
}

double getResid(double* x, double* b, const int size)
{
   int i;
   double localres,resmag;
   double global_resmag;
   
   // Prepare for async send/recv
   MPI_Request request[4];
   int requests;
   MPI_Status status[4];
   
   // grab the left and right buffer.
   double left_buffer = 0.0;
   double right_buffer = 0.0;
   
   requests=0;
      
   // Fill the left buffer. Send to the right, listen from the left.
   MPI_Isend(&x[size-1],   1, MPI_DOUBLE, (my_rank+1)%world_size, 1, MPI_COMM_WORLD, request + requests++);
   MPI_Irecv(&left_buffer, 1, MPI_DOUBLE, (my_rank+world_size-1)%world_size, 1, MPI_COMM_WORLD, request + requests++);

   //printf("I am rank %d of %d and I received %.8e from the left.\n", my_rank, world_size, left_buffer);

   // Fill the right buffer. Send to the left, listen from the right.
   MPI_Isend(&x[0],   1, MPI_DOUBLE, (my_rank+world_size-1)%world_size, 0, MPI_COMM_WORLD, request + requests++);
   MPI_Irecv(&right_buffer, 1, MPI_DOUBLE, (my_rank+1)%world_size, 0, MPI_COMM_WORLD, request + requests++);
   //printf("I am rank %d of %d and I received %.8e from the right.\n", my_rank, world_size, right_buffer);

   i = 0;
   localres = 0.0;
   global_resmag = 0.0;
   resmag = 0.0;

   // Loop over rest.
   for (i=1;i<size-1;i++)
   {
      localres = (b[i] - x[i] + 0.5*(x[i+1] + x[i-1]));
      localres = localres*localres;
      resmag = resmag + localres;
   }
   
   // Wait for async.
   MPI_Waitall ( requests, request, status );
   
   // Handle boundaries, acknowledging 0 bcs.
   // Impose zero bc.
   if (my_rank != 0)
   {
      localres = (b[0] - x[0] + 0.5*(x[1] + left_buffer));
      localres = localres*localres;
      resmag = resmag + localres;
   }
   //else
   //{
   //   localres = (b[0] - x[0] + 0.5*(x[1] + 0.0));
   //   localres = localres*localres;
   //   resmag = resmag + localres;
   //}

   // Impose zero bc.
   if (my_rank != world_size-1)
   {
      localres = (b[size-1] - x[size-1] + 0.5*(right_buffer + x[size-2]));
      localres = localres*localres;
      resmag = resmag + localres;
   }
   //else
   //{
   //   localres = (b[size-1] - x[size-1] + 0.5*(0.0 + x[size-2]));
   //   localres = localres*localres;
   //   resmag = resmag + localres;
   //}

   //printf("I am rank %d of %d and I have a local square residual %.8e.\n", my_rank, world_size, resmag);

   // Reduce. 
   MPI_Allreduce(&resmag, &global_resmag, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
   
   return sqrt(global_resmag);
}