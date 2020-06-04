#!/bin/bash
#PBS -l walltime=00:02:00,nodes=2:ppn=4
#PBS -N integral_job
#PBS -q batch
cd $PBS_O_WORKDIR
for i in {1..8}
do
	mpirun --hostfile $PBS_NODEFILE -np $i ./a.out 1000
	mpirun --hostfile $PBS_NODEFILE -np $i ./a.out 1000000
	mpirun --hostfile $PBS_NODEFILE -np $i ./a.out 100000000
done
