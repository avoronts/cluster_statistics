#!/bin/bash

#SBATCH -N 6
#SBATCH -t 2000
##SBATCH -p debug

mpirun -n 12 ./program.x 
#>& dump.out

