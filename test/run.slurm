#!/bin/bash

#SBATCH -N 1
#SBATCH -t 14-00:00:00

echo `which mpirun`
mpirun -np 12 ./program.x > /dev/null

