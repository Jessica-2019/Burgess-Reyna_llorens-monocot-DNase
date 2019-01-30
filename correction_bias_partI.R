### Scripts for DNAse-I seq correction bias. 
### Reyna-Llorens I, 2018


require(dplyr);
require(reshape);
require(reshape2);
library(GenomicRanges);
library(gtools);
library(mixtools);
require(gplots);

#Before running this code a SeqBias.txt file must be generated from the nkd DNAseI seq experiments. 

args <- commandArgs(trailingOnly = TRUE)
dgfs <- args[1]
chrom <- args[2]
genome <-args[3]
bam <-args[4]
prefix <-args[5]

##

dgf <-read.delim(dgfs, head=F, sep="\t")


dgf$size<-dgf$V3 - dgf$V2
dgf$sf <-as.factor(dgf$size)
head(dgf)
nrow(dgf)
lev<-unique(dgf$sf)
summary(lev)

### Adjust all bases to the nearest multiple of 6 and save as bed file
num = as.integer(dgf$size[1])
num
if((num %% 2) == 0) {
  dgf.12 <- dgf %>% filter(sf==10 |sf==12| sf==14)  %>% 
    mutate(V3= ifelse(sf==10, V3+2, V3),
           V3= ifelse(sf==14, V3-2, V3))
  
  dgf.18 <- dgf %>% filter(sf==16 |sf==18| sf==20)  %>% 
    mutate(V3= ifelse(sf==16, V3+2, V3),
           V3= ifelse(sf==20, V3-2, V3))
  
  dgf.24 <- dgf %>% filter(sf==22|sf==24)  %>% 
    mutate(V3= ifelse(sf==22, V3+2, V3))
  
  } else { dgf.12 <- dgf %>% filter(sf==11 | sf==13)  %>% 
    mutate(V3= ifelse(sf==11, V3+1, V3),
           V3= ifelse(sf==13, V3-1, V3))
  
  dgf.18 <- dgf %>% filter(sf==15 | sf==17|sf==19)  %>% 
    mutate(V3= ifelse(sf==15, V3+3, V3),
           V3= ifelse(sf==17, V3+1, V3),
           V3= ifelse(sf==19, V3-1, V3))
  
  dgf.24 <- dgf %>% filter(sf==21 | sf==23 | sf==25)  %>% 
    mutate(V3= ifelse(sf==21, V3+3, V3),
           V3= ifelse(sf==23, V3+1, V3),
           V3= ifelse(sf==25, V3-1, V3))}

head(dgf.12)

#write

write.table(dgf.12[,1:4], file=paste0(prefix,"dgf12.bed"), sep="\t",
            row.names=FALSE,quote = FALSE, col.names = FALSE) 
write.table(dgf.18[,1:4], file=paste0(prefix,"dgf18.bed"), sep="\t",
            row.names=FALSE,quote = FALSE, col.names = FALSE) 
write.table(dgf.24[,1:4], file=paste0(prefix,"dgf24.bed"), sep="\t",
            row.names=FALSE,quote = FALSE, col.names = FALSE) 

#extend 25 + 3 more and to a multiple of 6 to get the fasta file 
system(paste0("bedtools slop -i ", prefix,"dgf12.bed -g ", chrom, " -b 28 > ext.12.28.bed"))#28+12
system(paste0("bedtools slop -i ", prefix,"dgf18.bed -g ", chrom, " -b 28 > ext.18.28.bed"))#28+18  
system(paste0("bedtools slop -i ", prefix,"dgf24.bed -g ", chrom, " -b 28 > ext.24.28.bed")) #28+24


# Obtain fasta file from extended DGFs 
system(paste0("bedtools getfasta -fi ",genome, " -bed ext.12.28.bed -fo ext.12.28.fasta"))
system(paste0("bedtools getfasta -fi ",genome, " -bed ext.18.28.bed -fo ext.18.28.fasta"))
system(paste0("bedtools getfasta -fi ",genome, " -bed ext.24.28.bed -fo ext.24.28.fasta"))

# Created the coverage matrix using dnase_to_javatreeview (bed file should only have 3 columns)

system(paste0("awk -v OFS='\t' '{print $1, $2, $3}' ", prefix,"dgf12.bed > dgf12.1.bed"))
system(paste0("awk -v OFS='\t' '{print $1, $2, $3}' ", prefix,"dgf18.bed > dgf18.1.bed"))
system(paste0("awk -v OFS='\t' '{print $1, $2, $3}' ", prefix,"dgf24.bed > dgf24.1.bed"))

system(paste0("dnase_to_javatreeview.py -w 34  -o -a dgf18.1.bed ",bam, " 18.ext.cov.csv"))
system(paste0("dnase_to_javatreeview.py -w 34  -o -a dgf12.1.bed ",bam, " 12.ext.cov.csv"))
system(paste0("dnase_to_javatreeview.py -w 34  -o -a dgf24.1.bed ",bam, " 24.ext.cov.csv"))

##############################################################################################################
###Proceed with partII of the script