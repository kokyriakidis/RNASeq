# RNASeq Kallisto-Sleuth Pipeline for PBS script submission in Cluster

## DEPENDENCIES
```
# R packages that must be installed prior in the Cluster #
source("http://bioconductor.org/biocLite.R")
biocLite("rhdf5")
install.packages("devtools")
library("devtools")
devtools::install_github("pachterlab/sleuth")
```
## INSTALLATION 
```
git clone https://github.com/kokyriakidis/RNASeq.git
cd RNASeq
```
## PIPELINE - HOWTO
#### 1) Put your RNA-Seq fastq files in FASTQ_FILES folder. RENAME THEM SO AS THE PREFIX OF THE FASTQ FILE IS THE SAMPLE_ID'S NAME

```
eg. FASTQ_FILES/SRR1804790_1.fastq.gz if the sample_id name is SRR1804790
    FASTQ_FILES/SRR1804790_2.fastq.gz if the sample_id name is SRR1804790

    FASTQ_FILES/SRR1804791_1.fastq.gz if the sample_id name is SRR1804791
    FASTQ_FILES/SRR1804791_2.fastq.gz if the sample_id name is SRR1804791

    FASTQ_FILES/SRR1804792_1.fastq.gz if the sample_id name is SRR1804792
    FASTQ_FILES/SRR1804792_2.fastq.gz if the sample_id name is SRR1804792

    FASTQ_FILES/SRR1804793_1.fastq.gz if the sample_id name is SRR1804793
    FASTQ_FILES/SRR1804793_2.fastq.gz if the sample_id name is SRR1804793
    .....
```

#### 2) Create folders in Preproccesing folder FOR EACH SAMPLE. THE NAME OF EACH FOLDER MUST BE THE SAMPLE_ID's NAME

```
eg. Preproccesing/SRR1804790
    Preproccesing/SRR1804791
    Preproccesing/SRR1804792
    Preproccesing/SRR1804793
    ......
```    

#### 3) How to edit the metadata.txt file in the Metadata directory (MUST BE A TAB DELIMITED FILE)

```
#Imagine that we have 2 samples for each of 2 conditions (untreated and treated):

sample_id  | condition  | path
sample01  | untreated  | path1
sample02  | untreated  | path2
sample03  | treated  | path3
sample04  | treated  | path4

path1, path2 ... pathn is the kallisto's output files directory for each sample.

#Example metadata.txt file

sample_id	condition	path
SRR1804790	untreated   /mnt/scratchdir/home/kyriakidk/Preproccesing/SRR1804790/
SRR1804791	untreated   /mnt/scratchdir/home/kyriakidk/Preproccesing/SRR1804791/
SRR1804792	treated /mnt/scratchdir/home/kyriakidk/Preproccesing/SRR1804792/
SRR1804793	treated /mnt/scratchdir/home/kyriakidk/Preproccesing/SRR1804793/
```

#### 4) Edit RNASeq_PrePro.sh script's PBS parameters to your needs
```
#!/bin/bash
#PBS -N PrePro_R_${ERR}
#PBS -q see
#PBS -j oe
#PBS -M email
#PBS -m bea
#PBS -l nodes=1:ppn=40,walltime=10:00:00
```
#### 5) Edit RNASeq_Kallisto.sh script's PBS parameters to your needs
```
#!/bin/bash
#PBS -N RNASeq_Analysis
#PBS -q see
#PBS -j oe
#PBS -M email
#PBS -m bea
#PBS -l nodes=1:ppn=40,walltime=10:00:00
```
#### 6) Edit RNASeq_Sleuth.sh script's PBS parameters to your needs AND SPECIFY IN THE diff_exp_sleuth FUNCTION THE 2 CONDITIONS, EXACTLY AS STATED IN THE metadata.txt's CONDITION COLUMN
```
#!/bin/bash
#PBS -N RNASeq_Sleuth
#PBS -q see
#PBS -j oe
#PBS -M email
#PBS -m bea
#PBS -l nodes=1:ppn=40,walltime=10:00:00

diff_exp_sleuth <- function("condition1", "condition2") 

# eg. diff_exp_sleuth <- function("untreated","treated") 

```
#### 7) Run RNASeq_PrePro.sh in cluster FOR EACH sample_id
```
qsub -v ERR=sample_id RNASeq_PrePro.sh

eg. qsub -v ERR=SRR1804790 RNASeq_PrePro.sh
    qsub -v ERR=SRR1804791 RNASeq_PrePro.sh
    qsub -v ERR=SRR1804792 RNASeq_PrePro.sh
    qsub -v ERR=SRR1804793 RNASeq_PrePro.sh
    .......

```
#### 8) Run RNASeq_Index.sh in cluster
```
qsub -v ERR=sample_id RNASeq_Index.sh
```
#### 9) Run RNASeq_Kallisto.sh in cluster FOR EACH sample_id
```
qsub -v ERR=sample_id RNASeq_Kallisto.sh

eg. qsub -v ERR=SRR1804790 RNASeq_Kallisto.sh
    qsub -v ERR=SRR1804791 RNASeq_Kallisto.sh
    qsub -v ERR=SRR1804792 RNASeq_Kallisto.sh
    qsub -v ERR=SRR1804793 RNASeq_Kallisto.sh
    .......

```

#### 10) Run RNASeq_Sleuth.sh in cluster FOR EACH sample_id
```
qsub -v ERR=sample_id RNASeq_Sleuth.sh

eg. qsub -v ERR=SRR1804790 RNASeq_Sleuth.sh
    qsub -v ERR=SRR1804791 RNASeq_Sleuth.sh
    qsub -v ERR=SRR1804792 RNASeq_Sleuth.sh
    qsub -v ERR=SRR1804793 RNASeq_Sleuth.sh
    .......

```
#### 11) The analysis' results will be inside the Analysis directory 
