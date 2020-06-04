#!/bin/bash
#PBS -l walltime=00:10:00,nodes=1:ppn=4
#PBS -N delete_1
#PBS -q batch
cd $PBS_O_WORKDIR
./a.out 1
./a.out 2
./a.out 3
./a.out 4
./a.out 5
./a.out 6
