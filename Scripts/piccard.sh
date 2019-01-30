#!/bin/bash



genome_fa=$1
bam=$2
echo ${bam}

for f in ${bam};
do
	samtools sort ${bam} -@6 -o ${bam}.srt 
	samtools index ${bam}.srt


#java -jar picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=${bam} O=${bam}_no_dup M=${bam%.bam}_marked_dup_metrics.txt
#remove a minimal mapq score of <10 and plot the scores from the tagAlign file to assees whether more is needed

#samtools sort -@ 6 -n ${bam}_no_dup | samtools fixmate -r -O bam - - | samtools view -q 42 -F 1804 -@ 6 -O BAM -u - | samtools sort -o ${bam}_filt -
#samtools index ${bam}_filt
	java -jar ~/picard.jar EstimateLibraryComplexity I=${bam}.srt O=${bam}_lib_complex_metrics.txt &
    java -jar ~/picard.jar QualityScoreDistribution I=${bam}.srt O=${bam}_qual_score_dist.txt CHART=${bam}_qual_score_dist.pdf &	# this by default only looks at 500000 reads - so just plot from tagAlign file (faster)
	java -jar ~/picard.jar CollectGcBiasMetrics I=${bam}.srt O=${bam}_gc_bias_metrics.txt CHART=${bam_in%.bam}_gc_bias_metrics.pdf S=${bam_in%.bam}_summary_metrics.txt R=${genome_fa} &
	java -jar ~/picard.jar CollectInsertSizeMetrics I=${bam}.srt O=${bam}.window500.hist_data H=${bam}.window500.hist_graph.pdf W=500

done
