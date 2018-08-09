#!/bin/bash
#PBS -N RNASeq_Analysis
#PBS -q see
#PBS -j oe
#PBS -M email
#PBS -m bea
#PBS -l nodes=1:ppn=40,walltime=10:00:00

##Building an index
./kallisto/kallisto index -i ./Resources/hg38_index.idx ./Resources/Homo_sapiens.GRCh38.rna.fa.gz

##Quantification. As many times as the sample's number
./kallisto/kallisto quant -i ./Resources/hg38_index.idx -o ./Preproccesing/${ERR}/ -b 100 ./Preproccesing/${ERR}/${ERR}_1_qt.fastq.gz ./Preproccesing/${ERR}/${ERR}_2_qt.fastq.gz  #FOR SINGLE-END READS: --single -l 180 -s 20 reads_1.fastq.gz
