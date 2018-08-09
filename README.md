# RNASeq

## 1) How to edit the metadata.txt file in the Metadata directory (MUST BE A TAB DELIMITED FILE)

```
#Imagine that we have 3 samples for each of 2 conditions (untreated and treated):

sample_id  | condition  | path
sample01  | untreated  | path1
sample02  | untreated  | path2
sample03  | untreated  | path3
sample04  | treated  | path4
sample05  | treated  | path5
sample06  | treated  | path6

path1, path2 ... pathn is the path where the kallisto's output files are for each sample.

```
```
#Example metadata.txt file

sample_id	condition	path
SRR1804790	Control	/Users/kokyriakidis/Desktop/rnaseq/results/SRR1804790/kallisto/
SRR1804791	Control	/Users/kokyriakidis/Desktop/rnaseq/results/SRR1804791/kallisto/
SRR1804792	SOX15	/Users/kokyriakidis/Desktop/rnaseq/results/SRR1804792/kallisto/
SRR1804793	SOX15	/Users/kokyriakidis/Desktop/rnaseq/results/SRR1804793/kallisto/
```
