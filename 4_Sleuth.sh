#!/bin/bash
#PBS -N Sleuth
#PBS -q see
#PBS -j oe
#PBS -M email
#PBS -m bea
#PBS -l nodes=1:ppn=20,walltime=10:00:00

export PATH=$DENOVOFS/groups/denovo/anaconda2/bin:$PATH

cd $PBS_O_WORKDIR

source activate RNASeq

Rscript ARG_SCRIPT.R "Control" "SOX15" &> $PBS_JOBID.log

source deactivate
