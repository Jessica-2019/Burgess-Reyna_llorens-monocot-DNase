#Script for permutation tests and overlap analysis between Z.mays DHS, DGF and ChIPseq data
#Author: Reyna-Llorens I, 2018

library(regioneR)
library(GenomicRanges)
library(rtracklayer)

setwd("/media/iar28/7cb991c7-8087-4444-ad62-da3ae7c40dae/home/iar28/documents/Zmays_ChIPseq/dnase/dnase_chipseq_AGVP3/")
databed <- read.table("chip_all_no_scaffold.bed", sep = "\t",header = F)
dgfsbed <- read.table("zmwl_dgf_no_scaffold.bed", sep = "\t",header = F)


bedtest <- databed[sample(nrow(databed), 5000), ]
dgftest <- dgfsbed[sample(nrow(dgfsbed), 5000), ]


my.granges <- GRanges(seqnames = databed$V1, 
                      ranges = IRanges(start = databed$V2,
                      end = databed$V3,
                      names = databed$V5),
                       length = databed$V4)


my.dgf <- GRanges(seqnames = dgfsbed$V1, 
                      ranges = IRanges(start = dgfsbed$V2,
                      end = dgfsbed$V3))


head(my.granges)

require(dplyr)
zmchrom <-read.table("Chr.Zm.size", head=F)
nrow(zmchrom)

zmchrom <-zmchrom[order(zmchrom$V1),] 
zmchrom<-zmchrom[1:12,]
head(zmchrom)

zmchrom$start <-1
zmchrom<-zmchrom[,c(1,3,2)]
colnames(zmchrom)<-c("chr","start","end")
zmchrom$chr<-as.character(zmchrom$chr)

maize<- GRanges(seqnames=Rle(zmchrom$chr),
                ranges=IRanges(zmchrom$start,
                end=zmchrom$end),lenght=zmchrom$end)

seqlengths(maize) <-c("1"=301476924, "10"=149632204,
                      "2"=237917468,"3"=232245527,
                      "4"=242062272,"5"=217959525,
                      "6"=169407836,"7"=176826311,
                      "8"=175377492,"9"=157038028,
                      "Mt"=569630,"Pt"=140384)


genome(maize)<-"maize"
seqinfo(maize)
numOverlaps(my.granges, my.dgf, count.once=TRUE)
numOverlaps(randomizeRegions(my.granges,genome = maize), my.dgf, count.once=TRUE)

pt <- overlapPermTest(A=my.granges, B=my.dgf, ntimes=100, genome=maize)
plot(pt)
lz <- localZScore(pt=pt, A=A, B=B)
plot(lz)

