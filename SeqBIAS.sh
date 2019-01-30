#!/bin/bash
#calculate DNAse-I cutting bias for each species
# After runing this script please run "correction_bias_partI.R"


for i in *.bam;
do
 bamToBed -i ${i} > ${i::-12}.bed
 
 less ${i::-12}.bed |awk 'OFS="\t" {if($6=="+") print $1,$2-3,$2+3,$4,$5,$6}' > ${i::-12}.pos_strand_6bp.bed
 less ${i::-12}.bed |awk 'OFS="\t" {if($6=="-") print $1,$3-3,$3+3,$4,$5,$6}' > ${i::-12}.neg_strand_6bp.bed 
 cat ${i::-12}.pos_strand_6bp.bed ${i::-12}.neg_strand_6bp.bed > ${i::-12}_strand_6bp.bed 
 awk '($2 >0)' ${i::-12}_strand_6bp.bed > ${i::-12}strand_6bp.bed 
done
 

#Fasta from Bed 
echo "fastaFromBed Sb"
fastaFromBed -s -fi ../../Genomes/Sb/Sbicolor_255_v2.0.softmasked.fa -bed Sbstrand_6bp.bed -fo Sb_strand_6bp.bed.fa
echo "fastaFromBed Bd"
fastaFromBed -s -fi ../../Genomes/Bd/Bdistachyon_283_assembly_v2.0.fa -bed Bdstrand_6bp.bed -fo Bd_strand_6bp.bed.fa
echo "fastaFromBed Zm"
fastaFromBed -s -fi ../../Genomes/Zm/Zmays_284_AGPv3.softmasked.fa -bed Zmstrand_6bp.bed -fo Zm_strand_6bp.bed.fa


echo "frequencies Sb"
less Sb_strand_6bp.bed.fa | grep -v Chr | awk '{print toupper($0)}' | sort | uniq -c > Sb_6mers_counted.txt
echo "frequencies Bd"
less Bd_strand_6bp.bed.fa | grep -v Bd | awk '{print toupper($0)}' | sort | uniq -c > Bd_6mers_counted.txt
echo "frequencies Zm"
less Zm_strand_6bp.bed.fa | grep -v scaffold | awk '{print toupper($0)}' | sort | uniq -c > Zm_6mers_counted.txt
echo "the rest"

less Sb_6mers_counted.txt | awk '{print $2,$1,$1/1979472}' > Sb_6mers_counted_freq.txt
less Bd_6mers_counted.txt | awk '{print $2,$1,$1/3936233}' > Bd_6mers_counted_freq.txt
less Zm_6mers_counted.txt | awk '{print $2,$1,$1/1093752}' > Zm_6mers_counted_freq.txt


join Sb_6mers_counted_freq.txt Sb6mers.txt > Sb6mers_counted_freq.txt
less Sb6mers_counted_freq.txt |awk 'OFS="\t" {print $0,$3/$5}' > SbSeqBias.txt
awk -v OFS="\t" '$1=$1' SbSeqBias.txt| awk '{print $1, $6}' >Sb_SeqBias.txt

join Bd_6mers_counted_freq.txt Bd6mers.txt > Bd6mers_counted_freq.txt
less Bd6mers_counted_freq.txt |awk 'OFS="\t" {print $0,$3/$5}' > BdSeqBias.txt
awk -v OFS="\t" '$1=$1' BdSeqBias.txt| awk '{print $1, $6}' >Bd_SeqBias.txt

join Zm_6mers_counted_freq.txt Zm6mers.txt > Zm6mers_counted_freq.txt
less Zm6mers_counted_freq.txt |awk 'OFS="\t" {print $0,$3/$5}' > ZmSeqBias.txt
awk -v OFS="\t" '$1=$1' ZmSeqBias.txt| awk '{print $1, $6}' >Zm_SeqBias.txt



#######################################################################

## Rscript correction_bias_partI.R ../DGFs/DGF_Zm_BS.bed ../Files_Chr_Size/Chr.Zm.size ../../Genomes/Zm/Zmays_284_AGPv3.softmasked.fa ../Maize/Zm_WL.srt.bam Zm
## Rscript correction_bias_partI.R ../DGFs/DGF_Zm_WL.bed ../Files_Chr_Size/Chr.Zm.size ../../Genomes/Zm/Zmays_284_AGPv3.softmasked.fa ../Maize/Zm_BS.srt.bam Zm
