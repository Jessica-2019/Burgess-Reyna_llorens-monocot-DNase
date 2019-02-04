#!/bin/bash
# Need to use newer version of samtools (v1.2)
# Add to path
PATH=$PATH:/home/sjb287/bin/UCSC/:/home/sjb287/bin/MACS-master/bin/:/home/sjb287/bin/spp/R/:/home/sjb287/bin/
echo `export PATH`

# Get input options
S1_R1=$1
S1_R2=$2
S2_R1=$3
S2_R2=$4
S3_R1=$5
S3_R2=$6
BOWTIE=$7
GENOME=$8
NAME=$9
THREADS=10

#Define Names of Reps
N_Rep1="${NAME}_repf1_MQ42"
N_Rep2="${NAME}_repf2_MQ42"
N_Rep3="${NAME}_repf3_MQ42"

Rep_1="${NAME}.REP1"
Rep_2="${NAME}.REP2"
Rep_3="${NAME}.REP3"

#Define Names of Output Tag Files
TAGALIGN1="${N_Rep1}.mp.filt.srt.tagAlign.gz"
TAGALIGN2="${N_Rep2}.mp.filt.srt.tagAlign.gz"
TAGALIGN3="${N_Rep3}.mp.filt.srt.tagAlign.gz"

#Define name for optimal peak file
OPTIMAL="spp.optimal.${NAME}.POOLED.MACS_peaks.narrowPeak.regionPeak.gz"
SAMFILE="${NAME}.mp.sam"
FILT_BAM_PREFIX="${NAME}.mp.filt.srt"
PREFIX_TAGALIGN="${FILT_BAM_PREFIX}.tagAlign"
TAGALIGN="${PREFIX_TAGALIGN}.gz"
OUTPEAKS="${FILT_BAM_PREFIX}.MACS2_peaks.narrowPeak"
OUTBED="${OPTIMAL::-3}.bed"
T_BEDFILE="${FILT_BAM_PREFIX}.SPOT.bed"
SUMMARY="${NAME}.sum.qc"

##################################################
#Run QC Analysis on each rep
##################################################
echo "Step 1.1: Performing QC Analysis of Sample 1."
echo `. /home/sjb287/bin/Final_QC_analysis.sh ${THREADS} ${BOWTIE} ${S1_R1} ${S1_R2} ${GENOME} ${N_Rep1}`

echo "Step 1.2: Performing QC Analysis of Sample 2."
echo `. /home/sjb287/bin/Final_QC_analysis.sh  ${THREADS} ${BOWTIE} ${S2_R1} ${S2_R2} ${GENOME} ${N_Rep2} `

echo "Step 1.3: Performing QC Analysis of Sample 3."
echo `. /home/sjb287/bin/Final_QC_analysis.sh ${THREADS} ${BOWTIE} ${S3_R1} ${S3_R2} ${GENOME} ${N_Rep3}`

##################################################
#Run IDR analysis
##################################################
echo "Step 2: Performing Self IDR analysis."
echo `. /home/sjb287/bin/Final_IDR_analysis.sh ${TAGALIGN1} ${TAGALIGN2} ${TAGALIGN3} ${GENOME} ${NAME}`

#################################################
#Create Optimal Peak Set
#################################################
echo "Step 3: Creating Optimal Peak Set."
#Read in the IDR summary file
#Select the number of IDR peaks for each rep comparison and the pooled sample and sore in array 'ar'
#numPeaks_Rep1_and2 (column 1), numPeaks_Rep1_and3 (column 2), numPeaks_Rep2_and3 (column 3), numPeaks_POOLED (column 7)
ar=`sed -n '2p' ${NAME}.sum.IDR.txt | awk 'BEGIN{OFS="\n"}{print $1,$2,$3,$7}'`
IFS=$'\n'

#Calculate the optimal threshold
#optThresh=max(max_numPeaks_Rep,numPeaks_Rep0)
optThresh=`echo "${ar[*]}" | sort -nr | head -n1`

#Sort the peaks from the pooled sample p or qvalue and signal value, select the top number of peaks corresponding to the threshold
echo `zcat ${NAME}.POOLED.MACS_peaks.narrowPeak.regionPeak.gz | sort -k7nr,7nr | head -n ${optThresh} > ${OPTIMAL}`
echo `awk 'BEGIN{OFS="\t"};{print $1, $2, $3}' ${OPTIMAL} > ${OUTBED}`

##################################################
#Run QC Analysis Part II
##################################################
echo "Step 4: Performing QC analysis Part II."
echo `printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "NMAP" "PCB" "NSC" "RSC" "QTAG" "PEAKSFDR" "SPOT"> ${SUMMARY}`
echo `. /home/sjb287/bin/QC_analysis_part_II.sh ${N_Rep1} ${OUTBED} ${SUMMARY} {Rep_1}`
echo `. /home/sjb287/bin/QC_analysis_part_II.sh ${N_Rep2} ${OUTBED} ${SUMMARY} {Rep_2}`
echo `. /home/sjb287/bin/QC_analysis_part_II.sh ${N_Rep3} ${OUTBED} ${SUMMARY} {Rep_3}`

##################################################
#Run Saturation Analysis
##################################################
#echo "Step 5: Running Saturation Analysis."
#echo `. /home/sjb287/bin/Final_saturation_analysis.sh ${GENOME} ${NAME}`
#echo `. /home/sjb287/bin/Final_saturation_analysis_2.sh ${GENOME} ${NAME}`
