#!/bin/bash
#PBS -N Out_Align
#PBS -l walltime=48:00:00
#PBS -l nodes=1:ppn=PROCESSORS

cd $PBS_O_WORKDIR

python Output_SISRS.py WGS_Pairwise
