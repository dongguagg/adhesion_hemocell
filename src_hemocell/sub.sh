#!/bin/bash
#SBATCH -p amd_256
#SBATCH -N 1
#SBATCH -n 64
source /public23/soft/modules/module.sh
module load hdf5/1.10.4-parallel-intel17 cmake/3.17.0 gcc/7.3.0-wzm
export PATH=/public23/home/sca2326/codes/new-Hemocell/HemoCell-master/tools/packCells:$PATH
export PATH=/public23/home/sca2326/codes/new-Hemocell/HemoCell-master/examples/pipeflow:$PATH
mpirun -np 64 pipeflow config.xml
