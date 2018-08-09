#!/bin/bash
#PBS -N RNASeq_Sleuth
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

### DEA analysis with with sleuth using kallisto output in R ###

##Load function script inside R
source("./diff_exp_sleuth.R")

diff_exp_sleuth <- function(condition1, condition2)
