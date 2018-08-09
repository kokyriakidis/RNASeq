#!/bin/bash
#PBS -N PrePro_N_${ERR}
#PBS -q see
#PBS -j oe
#PBS -M email
#PBS -m bea
#PBS -l nodes=1:ppn=40,walltime=10:00:00

cd $PBS_O_WORKDIR

##Adapter + Filtering
./bbmap/bbduk.sh threads=$PBS_NP in1=./FASTQ_FILES/${ERR}_1.fastq.gz in2=./FASTQ_FILES/${ERR}_2.fastq.gz out1=./Preproccesing/${ERR}/${ERR}_1_trimmed.fastq.gz out2=./Preproccesing/${ERR}/${ERR}_2_trimmed.fastq.gz ref=./bbmap/resources/adapters.fa ktrim=r k=23 mink=10 hdist=1 tpe tbo ftm=5 ordered overwrite=true
mv ./Preproccesing/${ERR}/${ERR}_1_trimmed.fastq.gz ./Preproccesing/${ERR}/temp_1.fastq.gz
mv ./Preproccesing/${ERR}/${ERR}_2_trimmed.fastq.gz ./Preproccesing/${ERR}/temp_2.fastq.gz

##Contaminant Filtering
./bbmap/bbduk.sh threads=$PBS_NP in1=./Preproccesing/${ERR}/temp_1.fastq.gz in2=./Preproccesing/${ERR}/temp_2.fastq.gz out1=./Preproccesing/${ERR}/${ERR}_1_unmatched.fastq.gz out2=./Preproccesing/${ERR}/${ERR}_2_unmatched.fastq.gz outm1=./Preproccesing/${ERR}/${ERR}_1_matched.fastq.gz outm2=./Preproccesing/${ERR}/${ERR}_2_matched.fastq.gz ref=./bbmap/resources/phix174_ill.ref.fa.gz,./bbmap/resources/sequencing_artifacts.fa.gz,./bbmap/resources/pJET1.2.fa,./bbmap/resources/mtst.fa k=27 hdist=1 ordered cardinality overwrite=true
rm ./Preproccesing/${ERR}/temp_1.fastq.gz
rm ./Preproccesing/${ERR}/temp_2.fastq.gz
rm ./Preproccesing/${ERR}/${ERR}_1_matched.fastq.gz
rm ./Preproccesing/${ERR}/${ERR}_2_matched.fastq.gz
mv ./Preproccesing/${ERR}/${ERR}_1_unmatched.fastq.gz ./Preproccesing/${ERR}/temp_1.fastq.gz
mv ./Preproccesing/${ERR}/${ERR}_2_unmatched.fastq.gz ./Preproccesing/${ERR}/temp_2.fastq.gz


##Error-Correction-Phase-1
./bbmap/bbmerge.sh in1=./Preproccesing/${ERR}/temp_1.fastq.gz in2=./Preproccesing/${ERR}/temp_2.fastq.gz out1=./Preproccesing/${ERR}/${ERR}_erc1_1.fastq.gz out2=./Preproccesing/${ERR}/${ERR}_erc1_2.fastq.gz ecco mix vstrict ordered adapters=./bbmap/resources/adapters.fa overwrite=true
rm ./Preproccesing/${ERR}/temp_1.fastq.gz
rm ./Preproccesing/${ERR}/temp_2.fastq.gz
mv ./Preproccesing/${ERR}/${ERR}_erc1_1.fastq.gz ./Preproccesing/${ERR}/temp_1.fastq.gz
mv ./Preproccesing/${ERR}/${ERR}_erc1_2.fastq.gz ./Preproccesing/${ERR}/temp_2.fastq.gz


##Error-Correction-Phase-2
./bbmap/clumpify.sh in1=./Preproccesing/${ERR}/temp_1.fastq.gz in2=./Preproccesing/${ERR}/temp_2.fastq.gz out1=./Preproccesing/${ERR}/${ERR}_erc2_1.fastq.gz out2=./Preproccesing/${ERR}/${ERR}_erc2_2.fastq.gz ecc passes=4 reorder overwrite=true
rm ./Preproccesing/${ERR}/temp_1.fastq.gz
rm ./Preproccesing/${ERR}/temp_2.fastq.gz
mv ./Preproccesing/${ERR}/${ERR}_erc2_1.fastq.gz ./Preproccesing/${ERR}/temp_1.fastq.gz
mv ./Preproccesing/${ERR}/${ERR}_erc2_2.fastq.gz ./Preproccesing/${ERR}/temp_2.fastq.gz


##Error-Correction-Phase-3
./bbmap/tadpole.sh threads=$PBS_NP in1=./Preproccesing/${ERR}/temp_1.fastq.gz in2=./Preproccesing/${ERR}/temp_2.fastq.gz out1=./Preproccesing/${ERR}/${ERR}_erc3_1.fastq.gz out2=./Preproccesing/${ERR}/${ERR}_erc3_2.fastq.gz ecc k=62 ordered overwrite=true
rm ./Preproccesing/${ERR}/temp_1.fastq.gz
rm ./Preproccesing/${ERR}/temp_2.fastq.gz
mv ./Preproccesing/${ERR}/${ERR}_erc3_1.fastq.gz ./Preproccesing/${ERR}/temp_1.fastq.gz
mv ./Preproccesing/${ERR}/${ERR}_erc3_2.fastq.gz ./Preproccesing/${ERR}/temp_2.fastq.gz


##Quality-trim
./bbmap/bbduk.sh threads=$PBS_NP in1=./Preproccesing/${ERR}/temp_1.fastq.gz in2=./Preproccesing/${ERR}/temp_2.fastq.gz out1=./Preproccesing/${ERR}/${ERR}_1_qt.fastq.gz out2=./Preproccesing/${ERR}/${ERR}_2_qt.fastq.gz qtrim=r trimq=10 minlen=36 ordered overwrite=true
rm ./Preproccesing/${ERR}/temp_1.fastq.gz
rm ./Preproccesing/${ERR}/temp_2.fastq.gz
