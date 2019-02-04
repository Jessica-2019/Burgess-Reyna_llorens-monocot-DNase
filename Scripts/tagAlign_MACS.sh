#!/bin/bash
#Script to generate tagAlingn files and run MACS

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
##### second stage of bam processing to create tagAlign files #####
echo "tagAlign"
for f in *.srt.bam;
do
	if [ ! -f ${f%.srt.bam}.tagAlign ]; then
	samtools sort -@${THREADS} -n -o ${f}_n ${f};
	samtools view -@${THREADS} -b ${f}_n | bamToBed -i stdin -bedpe | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$6,"N",$8,$9}' > ${f%.srt.bam}.tagAlign;
	rm ${f}_n;
	fi;
done
echo "done tagAlign"
##### call macs2 to call DHS regions #####

echo "MACS"
#size of non- N or masked is 666303941 which is likely closer to the usable size (Ns=73515725)
SHIFT=-75
EXTSIZE=150
LLOCAL=50000

for f in *tagAlign;
do
	macs2 callpeak -t ${f} -f BED -g ${GENOME} -n ${f}.MACS -p 1e-1 --nomodel --extsize ${EXTSIZE} --shift ${SHIFT} --llocal ${LLOCAL};
done


