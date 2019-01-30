#!/bin/bash


# move bed files into one folder
for d in */*/; 
do
echo ${d}
samp=`basename ${d}`
cd ${d};
mv *.fastq.gz ../

cd ../..;

done

###########################

#cat files 

for d in Maize/;
do 
echo ${d}
samp=`basename ${d}`
cd ${d};
zcat *R1_001.fastq.gz |gzip -f > ${samp}_1.fastq.gz
zcat *R2_001.fastq.gz |gzip -f > ${samp}_2.fastq.gz
 
cd ../;

done

## Count the number of bases in a genome 

grep -v ">" file.fasta | wc | awk '{print $3-$1}'