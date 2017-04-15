#!/bin/sh
#
#$ -pe mpi_4_tasks_per_node 4
#
# Invoke mpirun.
# SGE sets $NSLOTS as the total number of processors
#
mpirun -np $NSLOTS ./my_2d_jacobi_mpi
