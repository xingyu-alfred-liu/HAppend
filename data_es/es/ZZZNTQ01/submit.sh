#!/bin/bash
#SBATCH -J relaxation # Job name
#SBATCH -n 56 # Number of total cores
#SBATCH -N 1 # Number of nodes
#SBATCH --mem-per-cpu=2048
#SBATCH -o j_%j.out # File to which STDOUT will be written %j is the job #
#SBATCH -p cpu
#SBATCH -t 06:00:00

echo "Job started on `hostname` at `date`"
mpirun -np 56 /home/xingyu/software/FHI-aims-05062019/aims.180607.scalapack.mpi.x &> aims.out
echo " "
echo "Job Ended at `date`"
