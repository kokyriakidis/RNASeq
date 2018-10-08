#!/bin/bash
#PBS -N Kallisto_${ERR}
#PBS -q see
#PBS -j oe
#PBS -M email
#PBS -m bea
#PBS -l nodes=1:ppn=20,walltime=10:00:00

cd $PBS_O_WORKDIR

##Quantification
./kallisto/kallisto quant -t $PBS_NP -i ./Resources/hg38_index.idx -o ./Preproccesing/${ERR}/ -b 100 ./Preproccesing/${ERR}/${ERR}_1_qt.fastq.gz ./Preproccesing/${ERR}/${ERR}_2_qt.fastq.gz

#FOR SINGLE-END READS: --single -l 180 -s 20 reads_1.fastq.gz
