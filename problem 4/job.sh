#!/bin/bash
#PBS -l walltime=00:10:00,nodes=1:ppn=4
#PBS -N WTF
#PBS -q batch
cd $PBS_O_WORKDIR
./a.out 2 0.0001 0.05
./a.out 4 0.0001 0.05
./a.out 6 0.0001 0.05
./a.out 8 0.0001 0.05
