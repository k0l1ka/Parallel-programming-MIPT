#!/bin/bash
#PBS -l walltime=00:10:00,nodes=3:ppn=4
#PBS -N cool_job
#PBS -q batch
cd $PBS_O_WORKDIR
for i in {1..12}
do
	mpirun -np $i --hostfile $PBS_NODEFILE ./a.out 2000 -1
done
for i in {1..12}
do
	mpirun -np $i --hostfile $PBS_NODEFILE ./a.out 10000 -1
done
for i in {1..12}
do
	mpirun -np $i --hostfile $PBS_NODEFILE ./a.out 50000 -1
done
