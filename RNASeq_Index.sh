#!/bin/bash
#PBS -N RNASeq_Index
#PBS -q see
#PBS -j oe
#PBS -M email
#PBS -m bea
#PBS -l nodes=1:ppn=40,walltime=10:00:00

cd $PBS_O_WORKDIR

##Building an index
./kallisto/kallisto index -i ./Resources/hg38_index.idx ./Resources/Homo_sapiens.GRCh38.rna.fa.gz
