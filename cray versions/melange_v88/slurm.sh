#!/bin/bash
#SBATCH -J ics432project # Name for your job
#SBATCH -n 1 # Number of tasks when using MPI. Default is 1
#SBATCH -c 20 # Number of cores requested, Default is 1 (total cores requested = tasks x cores)
#SBATCH -N 1 # Number of nodes to spread cores across - default is 1 - if you are not using MPI this should likely be 1
#SBATCH -t 1000 # Runtime in minutes
#SBATCH -p community.q# Partition to submit to the standard compute node partition in this example
#SBATCH -o results.out # Standard out goes to this file
#SBATCH -e results.out # Standard err goes to this file 
#SBATCH --mail-type ALL
#SBATCH --mail-user jeffreyt@hawaii.edu

for ((T=1; T <= 10; T++))
do
  echo "Test T=$T" >> results.out
  icc -DNUM_THREADS=20 -fopenmp -Ofast -c ezxml.c -o ezxml.o
  icc -DNUM_THREADS=20 -fopenmp -Ofast -c graphics.cpp -o graphics.o
  icc -DNUM_THREADS=20 -fopenmp -Ofast -c main.cpp -o main.o
  icc -DNUM_THREADS=20 -fopenmp -Ofast ezxml.o graphics.o main.o -o melange
  ./melange >> results.out
  echo "============================" >> results.out
done

