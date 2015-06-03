#!/bin/bash

#SBATCH -p quick
#SBATCH -N 1
#SBATCH -t 14-00:00:00

echo `which mpirun`
mpirun -np 6 ./program.x >& dump.out

