#!/bin/bash
#Script for aligning and cleanining DNAseI seq reads 
#Authors, Reyna-Llorens, Stevenson, Burgess

while getopts g:p:t:q: opt; do
	case $opt in
		g)
			if [ ${OPTARG::1} == "-" ]; then echo "No genome size given for -g. Genome size required for macs2 and should be the total genome length-N's"; exit 1; else GENOME=$OPTARG; fi
			;;
		t)
			echo "-t was triggered" >&2
			if [ ${OPTARG::1} == "-" ]; then echo "No thread number given for -t, using 8"; THREADS=8; else THREADS=$OPTARG; fi
			;;
		p)
			echo "-p was triggered" >&2
			if [ ${OPTARG::1} == "-" ]; then echo "No prefix given for -p, please ensure bowtie2-build has been used to build index"; exit 1; else PREFIX=$OPTARG; fi
			;;
		q)
			echo "-m was triggered" >&2
			if [ ${OPTARG::1} == "-" ]; then echo "No MAPQ threshold given for -q, using default of 30"; MAPQ=30; else MAPQ_THRESH=$OPTARG; fi	;;
	esac
done


##### mapping to create sam files #####
R1S=(*1.fastq*)
R2S=(*2.fastq*)
echo ${R1S[@]} 
echo ${R2S[@]}
x=0

BIF=`ls ${PREFIX}* | wc -l`
if [ ${BIF} -gt 2 ]; then echo "The genome files required for bowtie2 are here: `ls ${PREFIX}*`"
else echo "No bowtie2 index files found, please build using bowtie2-build"; exit 1; fi

echo "Runing bowtie2..."

for f in ${R1S[@]}; 
do
	if [ ! -f ${f%_R1*}.bam ]; then
	echo "running ${f} and ${R2S[x]}"
        
    /applications/bowtie2/bowtie2-2.2.6/bowtie2 -p $THREADS --local -D 15 -R 2 -N 0 -L 20 -i S,1,0.75 -q -x ${PREFIX} -1 ${f} -2 ${R2S[x]} | samtools view -@$THREADS -bS -|samtools  sort -@$THREADS -o ${f}.tmp2.bam
	samtools view -F 780 -h -f 2 -q 42 -@$THREADS -u ${f}.tmp2.bam | samtools sort -n -o ${f}.tmp3.bam -@$THREADS
	samtools fixmate -r -O bam ${f}.tmp3.bam ${f}.tmp4.bam 
	samtools sort -@$THREADS ${f}.tmp4.bam -o ${f}.tmp5.bam 
	rm ${f}.tmp2.bam
	rm ${f}.tmp3.bam
	rm ${f}.tmp4.bam
   	rm ${f}.tmp4.bam.bai
	samtools index ${f}.tmp5.bam
	samtools flagstat ${f}.tmp5.bam > ${f}.qc
	samtools sort -@$THREADS -n ${f}.tmp5.bam -o ${f}.bam_srt
	rm ${f}.tmp5.bam
	rm ${f}.tmp5.bam.bai
	x=$(expr $x + 1);
	fi;
done

echo "done filtering proceed with tagAlign"



