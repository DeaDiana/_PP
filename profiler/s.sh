#!/bin/sh
#PBS -l select=2:mpiprocs=2
cat $PBS_NODEFILES
cd $PBS_O_WORKDIR
mpirun -np 2 ./qw

//364749.hpcsuvir1.hpc.nusc.ru