#!/bin/bash

#SBATCH -N 6
#SBATCH -t 2000
##SBATCH -p debug

mpirun -n 6 ./program.x >& dump.out

