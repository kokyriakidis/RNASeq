# RNASeq Kallisto-Sleuth 

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

#### 3) How to edit the metadata.txt file in the Metadata directory (MUST BE A TAB DELIMITED FILE!)

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

#### 4) Edit 4_Sleuth.sh AND SPECIFY IN THE ARG_SCRIPT.R THE 2 CONDITIONS, EXACTLY AS STATED IN THE metadata.txt's CONDITION COLUMN
```
Rscript ARG_SCRIPT.R **_"Control" "SOX15"_** &> $PBS_JOBID.log 

DO NOT change "export PATH=$DENOVOFS/groups/denovo/anaconda2/bin:$PATH". You can ONLY change "cd $PBS_O_WORKDIR"

```

#### 5) Run 1_Index.sh in cluster
```
qsub 1_Index.sh
```

#### 6) Run 2_PrePro.sh in cluster FOR EACH sample_id
```
qsub -v ERR=sample_id 2_PrePro.sh

eg. qsub -v ERR=SRR1804790 2_PrePro.sh
    qsub -v ERR=SRR1804791 2_PrePro.sh
    qsub -v ERR=SRR1804792 2_PrePro.sh
    qsub -v ERR=SRR1804793 2_PrePro.sh
    .......

```

#### 7) Run 3_Kallisto.sh in cluster FOR EACH sample_id
```
qsub -v ERR=sample_id 3_Kallisto.sh

eg. qsub -v ERR=SRR1804790 3_Kallisto.sh
    qsub -v ERR=SRR1804791 3_Kallisto.sh
    qsub -v ERR=SRR1804792 3_Kallisto.sh
    qsub -v ERR=SRR1804793 3_Kallisto.sh
    .......

```

#### 8) Run 4_Sleuth.sh
```
qsub 4_Sleuth.sh

```
#### 9) The analysis' results will be inside the Analysis directory 
