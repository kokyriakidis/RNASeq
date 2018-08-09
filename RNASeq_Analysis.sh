#!/bin/bash
#PBS -N RNASeq_Analysis
#PBS -q see
#PBS -j oe
#PBS -M email
#PBS -m bea
#PBS -l nodes=1:ppn=40,walltime=10:00:00

module load R

###Libraries that must be installed in cluster###
#source("http://bioconductor.org/biocLite.R")
#biocLite("rhdf5")
#install.packages("devtools")
#library("devtools")
#devtools::install_github("pachterlab/sleuth")


##Building an index
./kallisto/kallisto index -i ./Resources/hg38_index.idx ./Resources/Homo_sapiens.GRCh38.rna.fa.gz

##Quantification. As many times as the sample's number
./kallisto/kallisto quant -i ./Resources/hg38_index.idx -o ./Preproccesing/${ERR1}/ -b 100 ./Preproccesing/${ERR1}/${ERR1}_1_qt.fastq.gz ./Preproccesing/${ERR1}/${ERR1}_2_qt.fastq.gz  #FOR SINGLE-END READS: --single -l 180 -s 20 reads_1.fastq.gz
./kallisto/kallisto quant -i ./Resources/hg38_index.idx -o ./Preproccesing/${ERR2}/ -b 100 ./Preproccesing/${ERR2}/${ERR2}_1_qt.fastq.gz ./Preproccesing/${ERR2}/${ERR2}_2_qt.fastq.gz  #FOR SINGLE-END READS: --single -l 180 -s 20 reads_1.fastq.gz
.....
./kallisto/kallisto quant -i ./Resources/hg38_index.idx -o ./Preproccesing/${ERRn}/ -b 100 ./Preproccesing/${ERRn}/${ERRn}_1_qt.fastq.gz ./Preproccesing/${ERRn}/${ERRn}_2_qt.fastq.gz  #FOR SINGLE-END READS: --single -l 180 -s 20 reads_1.fastq.gz


### DEA analysis with with sleuth using kallisto output in R ###

##Load function script inside R
source("./diff_exp_sleuth.R")

diff_exp_sleuth <- function(condition1, condition2)
