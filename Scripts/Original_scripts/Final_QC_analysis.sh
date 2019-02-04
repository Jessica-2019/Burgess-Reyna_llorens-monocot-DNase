#!/bin/bash
# Need to use newer version of samtools (v1.2)

# Add to path
PATH=$PATH:/home/sjb287/bin/UCSC/:/home/sjb287/bin/MACS-master/bin/:/home/sjb287/bin/spp/R/:/home/sjb287/bin/
echo `export PATH`

# Get input options
THREADS=$1
BOWTIE=$2 # Bowtie file
R1=$3
R2=$4
GENOME=$5
NAME=$6

# Define outfile
SAMFILE="${NAME}.mp.sam"
TEMPFILE1="${SAMFILE::-4}.bam"
TEMPFILE2="${SAMFILE::-4}.srt"
TEMPFILE3="$TEMPFILE2.bam"
SPP_OUT="${SAMFILE::-4}_SPP_STATS.tab"
FILT_BAM_PREFIX="${NAME}.mp.filt.srt"
FILT_BAM_PREFIX_n="${NAME}.mp.filt.nsrt"
FILT_BAM_FILE="${FILT_BAM_PREFIX}.bam"
FILT_BAM_FILE_n="${FILT_BAM_PREFIX_n}.bam"
FILT_BAM_INDEX_FILE="${FILT_BAM_PREFIX}.bai"
FILT_BAM_FILE_MAPSTATS="${FILT_BAM_PREFIX}.flagstat.qc" # QC file
TMP_FILT_BAM_PREFIX="tmp.${FILT_BAM_PREFIX}.nmsrt"
TMP_FILT_BAM_FILE="${TMP_FILT_BAM_PREFIX}.bam"
MAPQ_THRESH=42
PBC_FILE_QC="${FILT_BAM_PREFIX}.pbc.qc"
PREFIX_TAGALIGN="${FILT_BAM_PREFIX}.tagAlign"
TAGALIGN="${PREFIX_TAGALIGN}.gz"
TAGALIGNPDF="${PREFIX_TAGALIGN}.pdf"
PEAKFILE="${FILT_BAM_PREFIX}.MACS2"
OUTPEAKS="${FILT_BAM_PREFIX}.MACS2_peaks.narrowPeak"
OUTBED="${OUTPEAKS}.bed"
T_BEDFILE="${FILT_BAM_PREFIX}.SPOT.bed"
SUMMARY="${FILT_BAM_PREFIX}.sum.qc"


# Align reads
echo "Step 1: Aligning reads."
echo `bowtie2 -p $THREADS --local -D 15 -R 2 -N 0 -L 20 -i S,1,0.75 -q -x $BOWTIE -1 $R1 -2 $R2, -S $SAMFILE`

# Create sorted, indexed BAM file and clean up
#echo "Step 2: Creating sorted BAM file."
echo `samtools view -@11 -bS $SAMFILE > $TEMPFILE1`
echo `samtools sort -@11 $TEMPFILE1 $TEMPFILE2`

# =============================
# Remove unmapped, mate unmapped
# not primary alignment, reads failing platform
# Remove low MAPQ reads
# Only keep properly paired reads
# Obtain name sorted BAM file
# ==================
# See website for sam FLAG explanations: http://broadinstitute.github.io/picard/explain-flags.html
# Step 1: samtools view -F 780 -h -f 2 -q ${MAPQ_THRESH} -u ${RAW_BAM_FILE}
# -f = Only output reads with FLAG='2' which means (read mapped in proper pair)
# -F = Do not output reads with FLAG='780' which means (read unmapped, mate unmapped, not primary alignment, read fails platform QC)
# -q = Only include reads with a mapping quality better than threshold
# -h = output with BAM header
# -u = output is uncompressed BAM file which saves time when piping to stdout
# Step 2: samtools sort -n - ${TMP_FILT_BAM_PREFIX} # sort tempfile by name (-n) 
# Step 3: samtools fixmate -r -O bam ${TMP_FILT_BAM_FILE} # Fix mate pairs after MAPQ threshold check input needs to be name sorted
# -r = remove secondary or unmapped reads
#= -O = specifies output is BAM when piping to stdout
# Step 4: samtools view -F 1804 -f 2 -u 
# Step 5: samtools sort -n - ${FILT_BAM_PREFIX}
echo "Step 3: Filtering BAM file."
echo `samtools view -F 780 -h -f 2 -q ${MAPQ_THRESH} -@11 -u ${TEMPFILE3} | samtools sort -n - ${TMP_FILT_BAM_PREFIX} -@11`
echo `samtools fixmate -r -O bam ${TMP_FILT_BAM_FILE} - | samtools view -F 780 -f 2 -@11 -u - | samtools sort -@11 - ${FILT_BAM_PREFIX}`
echo `rm ${TMP_FILT_BAM_FILE}`

# Index Final BAM file needs position sorted reads
echo "Step 4: Indexing sorted, filtered BAM file."
echo `samtools index ${FILT_BAM_FILE}`
echo `samtools flagstat ${FILT_BAM_FILE} > ${FILT_BAM_FILE_MAPSTATS}`
echo `samtools sort -@12 -n ${FILT_BAM_FILE} ${FILT_BAM_PREFIX_n}`

##################################################
#Create tagAlignment file
##################################################
#bamToBed -bedpe form only reports one of two mate pair reads
echo `samtools view -@11 -b ${FILT_BAM_FILE_n} | bamToBed -i stdin -bedpe | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$6,"N",$8,$9}' | gzip -c > $TAGALIGN`
##################################################
#Calculate NSC and RSC
##################################################
echo "Step 6: Calculating NSC and RSC values."
#speak - it was necessary to input the predominant fragment length value.
#echo `Rscript /home/sjb287/bin/spp/R/run_spp.R -c=${TAGALIGN} -speak=74 -savp -p=10 -out=$SPP_OUT`
##################################################
# Calculate PCB - PCR bottleneck coefficient (although for DNase-Seq may have biological relevence rather than technical issue) 
##################################################
# Notes on PCB:
# A measure of library complexity, i.e. how skewed the distribution of read counts per location is towards 1 read per location.
# PBC = N1/Nd
# (where N1= number of genomic locations to which EXACTLY one unique mapping read maps, 
# and Nd = the number of genomic locations to which AT LEAST one unique mapping read maps, 
# i.e. the number of non-redundant, unique mapping reads).

# PBC is further described on the ENCODE Software Tools page. Provisionally, 0-0.5 is severe bottlenecking, 
# 0.5-0.8 is moderate bottlenecking, 0.8-0.9 is mild bottlenecking, while 0.9-1.0 is no bottlenecking. Very 
# low values can indicate a technical problem, such as PCR bias, or a biological finding, such as a very rare 
# genomic feature. Nuclease-based assays (DNase, MNase) detecting features with base-pair resolution (transcription 
# factor footprints, positioned nucleosomes) are expected to recover the same read multiple times, resulting in a 
# lower PBC score for these assays. Note that the most complex library, random DNA, would approach 1.0, thus the 
# very highest values can indicate technical problems with libraries. It is the practice for some labs outside of 
# ENCODE to remove redundant reads; after this has been done, the value for this metric is 1.0, and this metric 
# is not meaningful. 82% of TF ChIP, 89% of His ChIP, 77% of DNase, 98% of FAIRE, and 97% of control ENCODE datasets 
# have no or mild bottlenecking.
echo "Step 8: Calculating PCR bottleneck coefficient (PBC)."
echo `bedtools bamtobed -bedpe -i ${FILT_BAM_FILE_n} | awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$6,$9,$10}' | grep -v 'ChrM\ChrC'| sort | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' > ${PBC_FILE_QC}`
# Output columns
# TotalReadPairs 
# DistinctReadPairs 
# OneReadPair 
# TwoReadPairs 
# NRF=Distinct/Total 
# PBC1=OnePair/Distinct 
# PBC2=OnePair/TwoPair

#echo `rm ${SAMFILE}`
echo `rm ${TEMPFILE1}`
echo `rm ${TEMPFILE1}`
echo `rm ${TEMPFILE2}`
echo `rm ${TEMPFILE3}`
echo `rm ${FILT_BAM_FILE_n}`
