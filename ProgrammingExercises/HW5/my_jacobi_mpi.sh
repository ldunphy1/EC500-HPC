#!/bin/sh
#
#$ -pe mpi_4_tasks_per_node 4
#
# Invoke mpirun.
# SGE sets $NSLOTS as the total number of processors (8 for this example)
#
mpirun -np $NSLOTS ./my_jacobi_mpi