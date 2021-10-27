#!/bin/bash

#SBATCH --exclusive
#SBATCH -p RT_study
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 0-00:10 
#SBATCH --comment "ionwake"

./a.out > text4.txt
