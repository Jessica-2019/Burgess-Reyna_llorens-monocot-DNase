# Add to path
PATH=$PATH:/home/sjb287/bin/UCSC/:/home/sjb287/bin/MACS-master/bin/:/home/sjb287/bin/spp/R/:/home/sjb287/bin/
echo `export PATH`

#3 Sample inputs need to be in TAGALIGN format:
Rep1=$1 
Rep2=$2 
Rep3=$3
GENOME=$4
NAME=$5 

SUMMARY="${NAME}.sum.IDR.txt"
POOLED_Rep="${NAME}.POOLED.tagAlign.gz"
CHR="Chr.size"
GENOMESIZE=$(faSize ${GENOME} | head -1 | cut -d " " -f1)

#Names of pseduo reps
Rep1_pr1_ALIGN="${1::-12}.IDR.pr1.tagAlign.gz"
Rep1_pr2_ALIGN="${1::-12}.IDR.pr2.tagAlign.gz"
Rep2_pr1_ALIGN="${2::-12}.IDR.pr1.tagAlign.gz"
Rep2_pr2_ALIGN="${2::-12}.IDR.pr2.tagAlign.gz"
Rep3_pr1_ALIGN="${3::-12}.IDR.pr1.tagAlign.gz"
Rep3_pr2_ALIGN="${3::-12}.IDR.pr2.tagAlign.gz"
POOLED_pr1_ALIGN="${POOLED_Rep::-12}.IDR.pr1.tagAlign.gz"
POOLED_pr2_ALIGN="${POOLED_Rep::-12}.IDR.pr2.tagAlign.gz"

#MACS Output files for pooled samples
POOLED_MACS="${NAME}.POOLED.MACS"

#MACS Output files for individual samples
Rep1_MACS="${NAME}.REP1.MACS"
Rep2_MACS="${NAME}.REP2.MACS"
Rep3_MACS="${NAME}.REP3.MACS"
Rep1_MACSOUT="${Rep1_MACS}_peaks.narrowPeak"
Rep2_MACSOUT="${Rep2_MACS}_peaks.narrowPeak"
Rep3_MACSOUT="${Rep3_MACS}_peaks.narrowPeak"
POOLED_MACSOUT="${POOLED_MACS}_peaks.narrowPeak"

#MACS Output files for pseudo reps of individual samples
Rep1_pr1MACS="${NAME}.REP1.pr1.MACS"
Rep1_pr2MACS="${NAME}.REP1.pr2.MACS"
Rep1_pr1_MACSOUT="${Rep1_pr1MACS}_peaks.narrowPeak"
Rep1_pr2_MACSOUT="${Rep1_pr2MACS}_peaks.narrowPeak"

Rep2_pr1MACS="${NAME}.REP2.pr1.MACS"
Rep2_pr2MACS="${NAME}.REP2.pr2.MACS"
Rep2_pr1_MACSOUT="${Rep2_pr1MACS}_peaks.narrowPeak"
Rep2_pr2_MACSOUT="${Rep2_pr2MACS}_peaks.narrowPeak"

Rep3_pr1MACS="${NAME}.REP3.pr1.MACS"
Rep3_pr2MACS="${NAME}.REP3.pr2.MACS"
Rep3_pr1_MACSOUT="${Rep3_pr1MACS}_peaks.narrowPeak"
Rep3_pr2_MACSOUT="${Rep3_pr2MACS}_peaks.narrowPeak"

#MACS Output files for pseudo reps of the pooled sample
POOLED_pr1MACS="${NAME}.POOLED.pr1.MACS"
POOLED_pr2MACS="${NAME}.POOLED.pr2.MACS"
POOLED_pr1_MACSOUT="${POOLED_pr1MACS}_peaks.narrowPeak"
POOLED_pr2_MACSOUT="${POOLED_pr2MACS}_peaks.narrowPeak"

#Overlap file for
OVERLAP="${outputStub}r.overlap"
POOLED_OPTIMAL="${POOLED_MACSOUT}.final.optimal.bed"

#Step 1: Create pseduo replicates from individual samples
echo "Creating pseduo replicates"
echo `. /home/sjb287/bin/create_pseudo_rep.sh ${Rep1}`
echo `. /home/sjb287/bin/create_pseudo_rep.sh ${Rep2}`
echo `. /home/sjb287/bin/create_pseudo_rep.sh ${Rep3}`

#Step 2: Create pseduo replicates of pooled samples
echo "Pooling replicates"
echo `zcat ${REP1} ${Rep2} ${Rep3} | gzip -c > ${POOLED_Rep}`
echo "Creating pseduo replicates of pooled samples"
echo `. /home/sjb287/bin/create_pseudo_rep.sh ${POOLED_Rep}`

#Step 3: Call peaks on individual replicates
# You want to be calling >100K peaks if possible so the threshold is set very low -p 1e-1, this can be adjusted depending on read depth
echo "Calling peaks on individual samples."
echo `macs2 callpeak -t ${Rep1} -f BED -g ${GENOMESIZE} -n ${Rep1_MACS} -p 1e-1 --nomodel --extsize 150 --shift -75 --llocal 50000`
echo `macs2 callpeak -t ${Rep2} -f BED -g ${GENOMESIZE} -n ${Rep2_MACS} -p 1e-1 --nomodel --extsize 150 --shift -75 --llocal 50000`
echo `macs2 callpeak -t ${Rep3} -f BED -g ${GENOMESIZE} -n ${Rep3_MACS} -p 1e-1 --nomodel --extsize 150 --shift -75 --llocal 50000`

#Step 4: Call peaks on pooled replicates
echo "Calling peaks on pooled replicates"
echo `macs2 callpeak -t ${POOLED_Rep} -f BED -g ${GENOMESIZE} -n ${POOLED_MACS} -p 1e-1 --nomodel --extsize 150 --shift -75 --llocal 50000`

#Step 5: Call peaks on pseudo-reps from individual samples
echo "Calling peaks on pseudo-reps from individual samples"
echo `macs2 callpeak -t ${Rep1_pr1_ALIGN} -f BED -g ${GENOMESIZE} -n ${Rep1_pr1MACS} -p 1e-1 --nomodel --extsize 150 --shift -75 --llocal 50000`
echo `macs2 callpeak -t ${Rep1_pr2_ALIGN} -f BED -g ${GENOMESIZE} -n ${Rep1_pr2MACS} -p 1e-1 --nomodel --extsize 150 --shift -75 --llocal 50000`

echo `macs2 callpeak -t ${Rep2_pr1_ALIGN} -f BED -g ${GENOMESIZE} -n ${Rep2_pr1MACS} -p 1e-1 --nomodel --extsize 150 --shift -75 --llocal 50000`
echo `macs2 callpeak -t ${Rep2_pr2_ALIGN} -f BED -g ${GENOMESIZE} -n ${Rep2_pr2MACS} -p 1e-1 --nomodel --extsize 150 --shift -75 --llocal 50000`

echo `macs2 callpeak -t ${Rep3_pr1_ALIGN} -f BED -g ${GENOMESIZE} -n ${Rep3_pr1MACS} -p 1e-1 --nomodel --extsize 150 --shift -75 --llocal 50000`
echo `macs2 callpeak -t ${Rep3_pr2_ALIGN} -f BED -g ${GENOMESIZE} -n ${Rep3_pr2MACS} -p 1e-1 --nomodel --extsize 150 --shift -75 --llocal 50000`

#Step 6: Call peaks on pseudo-reps from pooled samples
echo "Calling peaks on pseduo reps from pooled samples"
echo `macs2 callpeak -t ${POOLED_pr1_ALIGN} -f BED -g ${GENOMESIZE} -n ${POOLED_pr1MACS} -p 1e-1 --nomodel --extsize 150 --shift -75 --llocal 50000`
echo `macs2 callpeak -t ${POOLED_pr2_ALIGN} -f BED -g ${GENOMESIZE} -n ${POOLED_pr2MACS} -p 1e-1 --nomodel --extsize 150 --shift -75 --llocal 50000`

#Step 7: 
# Run self-consistancy analysis
#Note from website:
#MACS2: MACS2 will output peaks in various formats (bed, xls and encodePeak). 
#You want to use the *encodePeak files. These are in the standard UCSC/ENCODE narrowPeak format. 
#Sort the *encodPeak files from best to worst using the -log10(pvalue) column i.e. column 8. 
#Then truncate the number of peaks to the top 100k-125k. Using more than this simply increases the running time of the IDR analysis with no advantage. 
#Infact using more peaks with MACS2 can cause problems with the IDR model because MACS2 seems to produce strange highly correlated peak scores for very 
#weak and noisy detections. This can confuse the IDR model.

echo "sorting and limiting number of called peaks to a max of 100K"
echo `sort -k 8nr,8nr ${Rep1_MACSOUT} | head -n 100000 | gzip -c > ${Rep1_MACSOUT}.regionPeak.gz`
echo `sort -k 8nr,8nr ${Rep2_MACSOUT} | head -n 100000 | gzip -c > ${Rep2_MACSOUT}.regionPeak.gz`
echo `sort -k 8nr,8nr ${Rep3_MACSOUT} | head -n 100000 | gzip -c > ${Rep3_MACSOUT}.regionPeak.gz`
echo `sort -k 8nr,8nr ${POOLED_MACSOUT} | head -n 100000 | gzip -c > ${POOLED_MACSOUT}.regionPeak.gz`

echo `sort -k 8nr,8nr ${Rep1_pr1_MACSOUT} | head -n 100000 | gzip -c > ${Rep1_pr1_MACSOUT}.regionPeak.gz`
echo `sort -k 8nr,8nr ${Rep1_pr2_MACSOUT} | head -n 100000 | gzip -c > ${Rep1_pr2_MACSOUT}.regionPeak.gz`

echo `sort -k 8nr,8nr ${Rep2_pr1_MACSOUT} | head -n 100000 | gzip -c > ${Rep2_pr1_MACSOUT}.regionPeak.gz`
echo `sort -k 8nr,8nr ${Rep2_pr2_MACSOUT} | head -n 100000 | gzip -c > ${Rep2_pr2_MACSOUT}.regionPeak.gz`

echo `sort -k 8nr,8nr ${Rep3_pr1_MACSOUT} | head -n 100000 | gzip -c > ${Rep3_pr1_MACSOUT}.regionPeak.gz`
echo `sort -k 8nr,8nr ${Rep3_pr2_MACSOUT} | head -n 100000 | gzip -c > ${Rep3_pr2_MACSOUT}.regionPeak.gz`

echo `sort -k 8nr,8nr ${POOLED_pr1_MACSOUT} | head -n 100000 | gzip -c > ${POOLED_pr1_MACSOUT}.regionPeak.gz`
echo `sort -k 8nr,8nr ${POOLED_pr2_MACSOUT} | head -n 100000 | gzip -c > ${POOLED_pr2_MACSOUT}.regionPeak.gz`

#Run Consistency Analysis Between Individual Samples
echo "Running IDR analysis"
echo `faSize $GENOME -detailed > ${CHR}`
echo "Running IDR analysis between samples"
echo `Rscript /home/sjb287/bin/spp/R/batch-consistency-analysis.r ${Rep1_MACSOUT}.regionPeak.gz ${Rep2_MACSOUT}.regionPeak.gz -1 0 F 0.02 ${CHR} "${NAME}.REP1vsREP2.r.output" "${NAME}.REP1vsREP2.r.overlap" "${NAME}.REP1vsREP2.npeaks.output" "${NAME}.REP1vsREP2.em.sav.output" "${NAME}.REP1vsREP2.uri.sav.output"`
echo `Rscript /home/sjb287/bin/spp/R/batch-consistency-analysis.r ${Rep1_MACSOUT}.regionPeak.gz ${Rep3_MACSOUT}.regionPeak.gz -1 0 F 0.02 ${CHR} "${NAME}.REP1vsREP3.r.output" "${NAME}.REP1vsREP3.r.overlap" "${NAME}.REP1vsREP3.npeaks.output" "${NAME}.REP1vsREP3.em.sav.output" "${NAME}.REP1vsREP3.uri.sav.output"`
echo `Rscript /home/sjb287/bin/spp/R/batch-consistency-analysis.r ${Rep2_MACSOUT}.regionPeak.gz ${Rep3_MACSOUT}.regionPeak.gz -1 0 F 0.02 ${CHR} "${NAME}.REP2vsREP3.r.output" "${NAME}.REP2vsREP3.r.overlap" "${NAME}.REP2vsREP3.npeaks.output" "${NAME}.REP2vsREP3.em.sav.output" "${NAME}.REP2vsREP3.uri.sav.output"`

#Run Consistency Analysis Between Pseudo Reps
echo "Running IDR analysis between pseudo reps"
echo `Rscript /home/sjb287/bin/spp/R/batch-consistency-analysis.r ${Rep1_pr1_MACSOUT}.regionPeak.gz ${Rep1_pr2_MACSOUT}.regionPeak.gz -1 0 F 0.02 ${CHR} "${NAME}.REP1.pr.r.output" "${NAME}.REP1.pr.r.overlap" "${NAME}.REP1.pr.npeaks.output" "${NAME}.REP1.pr.em.sav.output" "${NAME}.REP1.pr.uri.sav.output"`
echo `Rscript /home/sjb287/bin/spp/R/batch-consistency-analysis.r ${Rep2_pr1_MACSOUT}.regionPeak.gz ${Rep2_pr2_MACSOUT}.regionPeak.gz -1 0 F 0.02 ${CHR} "${NAME}.REP2.pr.r.output" "${NAME}.REP2.pr.r.overlap" "${NAME}.REP2.pr.npeaks.output" "${NAME}.REP2.pr.em.sav.output" "${NAME}.REP2.pr.uri.sav.output"`
echo `Rscript /home/sjb287/bin/spp/R/batch-consistency-analysis.r ${Rep3_pr1_MACSOUT}.regionPeak.gz ${Rep3_pr2_MACSOUT}.regionPeak.gz -1 0 F 0.02 ${CHR} "${NAME}.REP3.pr.r.output" "${NAME}.REP3.pr.r.overlap" "${NAME}.REP3.pr.npeaks.output" "${NAME}.REP3.pr.em.sav.output" "${NAME}.REP3.pr.uri.sav.output"`
echo `Rscript /home/sjb287/bin/spp/R/batch-consistency-analysis.r ${POOLED_pr1_MACSOUT}.regionPeak.gz ${POOLED_pr2_MACSOUT}.regionPeak.gz -1 0 F 0.02 ${CHR} "${NAME}.POOLED.pr.r.output" "${NAME}.POOLED.pr.r.overlap" "${NAME}.POOLED.pr.npeaks.output" "${NAME}.POOLED.pr.em.sav.output" "${NAME}.POOLED.pr.uri.sav.output"`


#For pooled-consistency analysis
#If you started with ~150 to 300K relaxed pre-IDR peaks for large genomes (human/mouse), then threshold of 0.0025 or 0.005 generally works well. 
#We use a tighter threshold for pooled-consistency since pooling and subsampling equalizes the pseudo-replicates in terms of data quality. 
#So we err on the side of caution and use more stringent thresholds. The equivalence between a pooled-consistency threshold of 0.0025 and 
#original replicate consistency threshold of 0.01 was calibrated based on a gold-standard pair of high quality replicate datasets for the CTCF transcription factor in human.
#You can also use 0.01 if you don't want to be very stringent.
#If you started with < 100K pre-IDR peaks for large genomes (human/mouse) , then a threshold of 0.01 is more appropriate. This is because the IDR sees a smaller 
#noise component and the IDR scores get weaker so we have to relax the thresholds. This is typically for use with peak callers that are unable to be adjusted to call large number of peaks (eg. PeakSeq or QuEST)

#For smaller genomes such as worm, if you start with ~15K to 40K peaks then once again IDR thresholds of 0.01 work well. 
#Each of these comparisons will give us a certain number of peaks that pass their respective IDR thresholds. We will refer to them as follows

#Original Replicate Threshold
numPeaks_Rep1_Rep2=$( awk '$11 <= 0.01 {print $0}' ${NAME}.REP1vsREP2.r.overlap | wc -l )
numPeaks_Rep1_Rep3=$( awk '$11 <= 0.01 {print $0}' ${NAME}.REP1vsREP3.r.overlap | wc -l )
numPeaks_Rep2_Rep3=$( awk '$11 <= 0.01 {print $0}' ${NAME}.REP2vsREP3.r.overlap | wc -l )

#SELF-CONSISTENCY THRESHOLDS
numPeaks_Rep1_pr=$( awk '$11 <= 0.01 {print $0}' ${NAME}.REP1.pr.r.overlap | wc -l )
numPeaks_Rep2_pr=$( awk '$11 <= 0.01 {print $0}' ${NAME}.REP2.pr.r.overlap | wc -l )
numPeaks_Rep3_pr=$( awk '$11 <= 0.01 {print $0}' ${NAME}.REP3.pr.r.overlap | wc -l )

#POOLED-PSEUDOREPLICATE THRESHOLD
numPeaks_Rep0=$( awk '$11 <= 0.0025 {print $0}' ${NAME}.POOLED.pr.r.overlap | wc -l )

echo `printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "numPeaks_Rep1_Rep2" "numPeaks_Rep1_Rep3" "numPeaks_Rep2_Rep3" "numPeaks_Rep1_pr" "numPeaks_Rep2_pr" "numPeaks_Rep3_pr" "numPeaks_Rep0" > ${SUMMARY}`
echo `printf "%s\t%s\t%s\t%s\t%s\t%s\t%s" "$numPeaks_Rep1_Rep2" "$numPeaks_Rep1_Rep3" "$numPeaks_Rep2_Rep3" "$numPeaks_Rep1_pr" "$numPeaks_Rep2_pr" "$numPeaks_Rep3_pr" "$numPeaks_Rep0" >> ${SUMMARY}`
