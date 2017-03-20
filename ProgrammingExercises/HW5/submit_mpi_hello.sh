#!/bin/sh
#
#$ -pe mpi_4_tasks_per_node 8
#
# Invoke mpirun.
# SGE sets $NSLOTS as the total number of processors (8 for this example)
#
mpirun -np $NSLOTS ./mpi_hello