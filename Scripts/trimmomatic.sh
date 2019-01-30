#!/bin/bash

#run trimmomatic on each sample
R1S=(*1.fastq*)
R2S=(*2.fastq*)
echo ${R1S[@]} 
echo ${R2S[@]}
x=0

for f in ${R1S[@]}; 
do

echo "Trimming sample ${d}"
java -jar ~/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 8 -phred33 ${f} ${R2S[x]} ${f::-16}.1.fq.gz ${f::-16}.1.unpaired.fq.gz ${f::-16}.2.fq.gz ${f::-16}.2.unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50
x=$(expr $x + 1);
done

