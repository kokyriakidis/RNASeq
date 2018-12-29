#!/bin/bash
#PBS -N PrePro_${ERR}
#PBS -q see
#PBS -j oe
#PBS -M email
#PBS -m bea
#PBS -l nodes=1:ppn=20,walltime=10:00:00

cd $PBS_O_WORKDIR

##Trim adapters.
./bbmap/bbduk.sh threads=$PBS_NP \
in1=./FASTQ_FILES/${ERR}_1.fastq.gz in2=./FASTQ_FILES/${ERR}_2.fastq.gz \
out1=./Preproccesing/${ERR}/${ERR}_1_trimmed.fastq.gz out2=./Preproccesing/${ERR}/${ERR}_2_trimmed.fastq.gz \
ref=./bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo ftm=5 ordered overwrite=true
mv ./Preproccesing/${ERR}/${ERR}_1_trimmed.fastq.gz ./Preproccesing/${ERR}/temp_1.fastq.gz
mv ./Preproccesing/${ERR}/${ERR}_2_trimmed.fastq.gz ./Preproccesing/${ERR}/temp_2.fastq.gz

##Remove synthetic artifacts and spike-ins by kmer-matching.
./bbmap/bbduk.sh threads=$PBS_NP \
in1=./Preproccesing/${ERR}/temp_1.fastq.gz in2=./Preproccesing/${ERR}/temp_2.fastq.gz \
out1=./Preproccesing/${ERR}/${ERR}_1_unmatched.fastq.gz out2=./Preproccesing/${ERR}/${ERR}_2_unmatched.fastq.gz \
ref=./bbmap/resources/phix174_ill.ref.fa.gz,./bbmap/resources/sequencing_artifacts.fa.gz,./bbmap/resources/pJET1.2.fa,./bbmap/resources/mtst.fa \
k=28 hdist=1 ordered cardinality overwrite=true
rm ./Preproccesing/${ERR}/temp_1.fastq.gz
rm ./Preproccesing/${ERR}/temp_2.fastq.gz
mv ./Preproccesing/${ERR}/${ERR}_1_unmatched.fastq.gz ./Preproccesing/${ERR}/temp_1.fastq.gz
mv ./Preproccesing/${ERR}/${ERR}_2_unmatched.fastq.gz ./Preproccesing/${ERR}/temp_2.fastq.gz

##Error Correction + Quality Correction.
fastp --thread=$PBS_NP \
--in1 ./Preproccesing/${ERR}/temp_1.fastq.gz --in2 ./Preproccesing/${ERR}/temp_2.fastq.gz \
--out1 ./Preproccesing/${ERR}/${ERR}_erc_1.fastq.gz --out2 ./Preproccesing/${ERR}/${ERR}_erc_2.fastq.gz \
--disable_adapter_trimming --disable_trim_poly_g \
--correction \
--disable_length_filtering \
--qualified_quality_phred 15 --unqualified_percent_limit 5 --n_base_limit 5
rm ./Preproccesing/${ERR}/temp_1.fastq.gz
rm ./Preproccesing/${ERR}/temp_2.fastq.gz
mv ./Preproccesing/${ERR}/${ERR}_erc_1.fastq.gz ./Preproccesing/${ERR}/temp_1.fastq.gz
mv ./Preproccesing/${ERR}/${ERR}_erc_2.fastq.gz ./Preproccesing/${ERR}/temp_2.fastq.gz

#Quality-trim + Length Filtering
./bbmap/bbduk.sh threads=$PBS_NP \
in1=./Preproccesing/${ERR}/temp_1.fastq.gz in2=./Preproccesing/${ERR}/temp_2.fastq.gz \
out1=./Preproccesing/${ERR}/${ERR}_qt_1.fastq.gz out2=./Preproccesing/${ERR}/${ERR}_qt_2.fastq.gz \
qtrim=rl trimq=15 minlen=45 ordered overwrite=true
rm ./Preproccesing/${ERR}/temp_1.fastq.gz
rm ./Preproccesing/${ERR}/temp_2.fastq.gz

source deactivate
