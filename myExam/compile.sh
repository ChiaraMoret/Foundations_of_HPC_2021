#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --nodes=1                
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:05:00
#SBATCH --job-name=compile
#SBATCH --error=compile.err
#SBATCH --output=complile.txt

#the aim of this bash script is to compile the programs on ORFEO
#therefore we need only very limited resources

#compile MPI
module load openMPI/4.1.5/gnu/12.2.1
mpicc MPIdensity.c -o MPIdensity.x -lm

#compile openMP
gcc -fopenmp oMPdensity.c  -o oMPdensity.x -lm


#compile auxiliary file
gcc readBinary.c -o readBinary.x -lm
gcc particles.c -o particles.x -lm
