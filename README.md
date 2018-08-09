# RNASeq

#### 1) Put your RNASeq fastq files in FASTQ_FILES folder. RENAME THEM SO AS THE PREFIX OF THE FASTQ FILE IS THE SAMPLE_ID'S NAME

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
SRR1804790	Control	/mnt/scratchdir/home/kyriakidk/Preproccesing/SRR1804790/
SRR1804791	Control	/mnt/scratchdir/home/kyriakidk/Preproccesing/SRR1804791/
SRR1804792	SOX15	/mnt/scratchdir/home/kyriakidk/Preproccesing/SRR1804792/
SRR1804793	SOX15	/mnt/scratchdir/home/kyriakidk/Preproccesing/SRR1804793/
```
