#!/bin/bash
# Add to path
PATH=$PATH:/home/sjb287/bin/UCSC/:/home/sjb287/bin/MACS-master/bin/:/home/sjb287/bin/spp/R/:/home/sjb287/bin/
echo `export PATH`

GENOME=$1
NAME=$2
TAGALIGN="${NAME}.POOLED.tagAlign.gz"
CHR="Chr.size"
GENOMESIZE=$(faSize ${GENOME} | head -1 | cut -d " " -f1)
OUTFILE="${NAME}.saturation.IDR.qc"
CHR="Chr.size"
GENOMESIZE=$(faSize ${GENOME} | head -1 | cut -d " " -f1)

#Calculate the total number of reads
TOT_READ=`zcat ${TAGALIGN}| wc -l | cut -d" " -f1`
#
#Create subsets of the reads
#run a for loop in bash, setting a counter to the value of anything between 0.05 and 1, at a set size of 0.05
#determine the number of reads required (V) by multiplying the total number of peaks by the counter value and rounding 
#shuffle all the reads and take a subset of the correct size 

for i in `seq 0.02 0.02 0.1`; do 

#Name output files
SUBSAMPLE="${NAME}.SUBSAMPLE.${i}"
SUBSAMPLE_MACS="${SUBSAMPLE}.MACS"
SUBSAMPLE_MACSOUT="${SUBSAMPLE_MACS}_peaks.narrowPeak"
SUBSAMPLE_pr1_MACS="${SUBSAMPLE}.pr1.MACS"
SUBSAMPLE_pr2_MACS="${SUBSAMPLE}.pr2.MACS"
SUBSAMPLE_pr1_MACSOUT="${SUBSAMPLE_pr1_MACS}_peaks.narrowPeak"
SUBSAMPLE_pr2_MACSOUT="${SUBSAMPLE_pr2_MACS}_peaks.narrowPeak"
SUBSAMPLE_pr1_ALIGN="${SUBSAMPLE}.IDR.pr1.tagAlign.gz"
SUBSAMPLE_pr2_ALIGN="${SUBSAMPLE}.IDR.pr2.tagAlign.gz"

#Create subsample
#Calculate number of reads
V=`awk 'BEGIN{printf "%.0f", '${TOT_READ}' * '${i}'}'`
echo `zcat ${TAGALIGN} | shuf | head -n $V | gzip ->${SUBSAMPLE}.tagAlign.gz`

#Create pseudoreps
echo `. /home/sjb287/bin/create_pseudo_rep.sh ${SUBSAMPLE}.tagAlign.gz`

#Run MACS2 on pseudoreps
echo `macs2 callpeak -t ${SUBSAMPLE_pr1_ALIGN} -f BED -g ${GENOMESIZE} -n ${SUBSAMPLE_pr1_MACS} -p 1e-1 --nomodel --extsize 150 --shift -75 --llocal 50000`
echo `macs2 callpeak -t ${SUBSAMPLE_pr2_ALIGN} -f BED -g ${GENOMESIZE} -n ${SUBSAMPLE_pr2_MACS} -p 1e-1 --nomodel --extsize 150 --shift -75 --llocal 50000`

#Sort peaks and take only top 100k
echo `sort -k 8nr,8nr ${SUBSAMPLE_pr1_MACSOUT} | head -n 100000 | gzip -c > ${SUBSAMPLE_pr1_MACSOUT}.regionPeak.gz`
echo `sort -k 8nr,8nr ${SUBSAMPLE_pr2_MACSOUT} | head -n 100000 | gzip -c > ${SUBSAMPLE_pr2_MACSOUT}.regionPeak.gz`

#Run batch consistancy analysis
echo `Rscript /home/sjb287/bin/spp/R/batch-consistency-analysis.r ${SUBSAMPLE_pr1_MACSOUT}.regionPeak.gz ${SUBSAMPLE_pr2_MACSOUT}.regionPeak.gz -1 0 F 0.02 ${CHR} "${SUBSAMPLE}.pr.r.output" "${SUBSAMPLE}.pr.r.overlap" "${SUBSAMPLE}.pr.npeaks.output" "${SUBSAMPLE}.pr.em.sav.output" "${SUBSAMPLE}.pr.uri.sav.output"`

numPeaks_Rep0=$( awk '$11 <= 0.0025 {print $0}' ${SUBSAMPLE}.pr.r.overlap | wc -l )
echo `printf "%s\t%s\n" "$V" "$numPeaks_Rep0" >> ${OUTFILE}`

done
for i in `seq 0.1 0.05 1`; do 
#Name output files
SUBSAMPLE="${NAME}.SUBSAMPLE.${i}"
SUBSAMPLE_MACS="${SUBSAMPLE}.MACS"
SUBSAMPLE_MACSOUT="${SUBSAMPLE_MACS}_peaks.narrowPeak"
SUBSAMPLE_pr1_MACS="${SUBSAMPLE}.pr1.MACS"
SUBSAMPLE_pr2_MACS="${SUBSAMPLE}.pr2.MACS"
SUBSAMPLE_pr1_MACSOUT="${SUBSAMPLE_pr1_MACS}_peaks.narrowPeak"
SUBSAMPLE_pr2_MACSOUT="${SUBSAMPLE_pr2_MACS}_peaks.narrowPeak"
SUBSAMPLE_pr1_ALIGN="${SUBSAMPLE}.IDR.pr1.tagAlign.gz"
SUBSAMPLE_pr2_ALIGN="${SUBSAMPLE}.IDR.pr2.tagAlign.gz"

#Create subsample
#Calculate number of reads
V=`awk 'BEGIN{printf "%.0f", '${TOT_READ}' * '${i}'}'`
echo `zcat ${TAGALIGN} | shuf | head -n $V | gzip ->${SUBSAMPLE}.tagAlign.gz`

#Create pseudoreps
echo `. /home/sjb287/bin/create_pseudo_rep.sh ${SUBSAMPLE}.tagAlign.gz`

#Run MACS2 on pseudoreps
echo `macs2 callpeak -t ${SUBSAMPLE_pr1_ALIGN} -f BED -g ${GENOMESIZE} -n ${SUBSAMPLE_pr1_MACS} -p 1e-1 --nomodel --extsize 150 --shift -75 --llocal 50000`
echo `macs2 callpeak -t ${SUBSAMPLE_pr2_ALIGN} -f BED -g ${GENOMESIZE} -n ${SUBSAMPLE_pr2_MACS} -p 1e-1 --nomodel --extsize 150 --shift -75 --llocal 50000`

#Sort peaks and take only top 100k
echo `sort -k 8nr,8nr ${SUBSAMPLE_pr1_MACSOUT} | head -n 100000 | gzip -c > ${SUBSAMPLE_pr1_MACSOUT}.regionPeak.gz`
echo `sort -k 8nr,8nr ${SUBSAMPLE_pr2_MACSOUT} | head -n 100000 | gzip -c > ${SUBSAMPLE_pr2_MACSOUT}.regionPeak.gz`

#Run batch consistancy analysis
echo `Rscript /home/sjb287/bin/spp/R/batch-consistency-analysis.r ${SUBSAMPLE_pr1_MACSOUT}.regionPeak.gz ${SUBSAMPLE_pr2_MACSOUT}.regionPeak.gz -1 0 F 0.02 ${CHR} "${SUBSAMPLE}.pr.r.output" "${SUBSAMPLE}.pr.r.overlap" "${SUBSAMPLE}.pr.npeaks.output" "${SUBSAMPLE}.pr.em.sav.output" "${SUBSAMPLE}.pr.uri.sav.output"`

numPeaks_Rep0=$( awk '$11 <= 0.0025 {print $0}' ${SUBSAMPLE}.pr.r.overlap | wc -l )
echo `printf "%s\t%s\n" "$V" "$numPeaks_Rep0" >> ${OUTFILE}`

done
