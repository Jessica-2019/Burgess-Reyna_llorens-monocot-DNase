##CORRECTION bias part II 
#runs MixtureModel.R
#copy MixtureModel.R, RebuildSignal.pl and SeqBias files in the folder
#modify MixtureModel.R so that the corresponding SeqBias file is used. 

require(dplyr);
require(reshape);
require(reshape2);
library(GenomicRanges);
library(gtools);
library(mixtools);
source("../../MixtureModel.r");

args <- commandArgs(trailingOnly = TRUE)
dhss <- args[1]
dgfs <-args[2]
prefix <-args[3]


f12.28 <-"ext.12.28.fasta"
f18.28 <-"ext.18.28.fasta"
f24.28 <-"ext.24.28.fasta"
dhs <-read.delim(dhss, head=F, sep="\t")
cov12 <- read.csv("12.ext.cov.csv", head=T, sep ="\t")
cov12<-cov12[,-c(1:3)]
cov18 <- read.csv("18.ext.cov.csv", head=T, sep ="\t")
cov18<-cov18[,-c(1:3)]
cov24 <- read.csv("24.ext.cov.csv", head=T, sep ="\t")
cov24<-cov24[,-c(1:3)]

dgf <-read.delim(dgfs, head=F, sep="\t")
dgf$size<-dgf$V3 - dgf$V2
dgf$sf <-as.factor(dgf$size)

### Adjust all bases to the nearest multiple of 6 and save as bed file
dgf.12 <- dgf %>% filter(sf==11 | sf==13)  %>% 
  mutate(V3= ifelse(sf==11, V3+1, V3),
         V3= ifelse(sf==13, V3-1, V3))

dgf.18 <- dgf %>% filter(sf==15 | sf==17|sf==19)  %>% 
  mutate(V3= ifelse(sf==15, V3+3, V3),
         V3= ifelse(sf==17, V3+1, V3),
         V3= ifelse(sf==19, V3-1, V3))

dgf.24 <- dgf %>% filter(sf==21 | sf==23 | sf==25)  %>% 
  mutate(V3= ifelse(sf==21, V3+3, V3),
         V3= ifelse(sf==23, V3+1, V3),
         V3= ifelse(sf==25, V3-1, V3))


m.12 <- MultMMixture_Full(TF_Bed=dgf.12[,1:3], Cuts=cov12, peakbed=dhs, Plot=T,
                          PadLen=22,Collapse=T,k=2,ReturnPar=T,Fixed=T,
                          Background ="Seq",FastaName=f12.28)

summary(m.12$llr)

m.24 <- MultMMixture_Full(TF_Bed=dgf.24[,1:3], Cuts=cov24, peakbed=dhs, Plot=T,
                          PadLen=28,Collapse=T,k=2,ReturnPar=T,Fixed=T,
                          Background ="Seq",FastaName=f24.28)


m.18 <- MultMMixture_Full(TF_Bed=dgf.18[,1:3], Cuts=cov18, peakbed=dhs, Plot=T,
                          PadLen=25,Collapse=T,k=2,ReturnPar=T,Fixed=T,
                          Background ="Seq",FastaName=f18.28)

summary(m.18$llr)


summary(m.24$llr)

#merge dgfnkd files with their corresponding Likelihoods and save it as a file. 

dgf.12$Likelihood <-m.12$llr
dgf.18$Likelihood <-m.18$llr
dgf.24$Likelihood <-m.24$llr
normdgf <-rbind(dgf.12, dgf.18, dgf.24)
normdgf <-normdgf[,-c(7,8)]
normdgf %>% filter (Likelihood > 0) %>% nrow() 

write.table(normdgf, sep="\t", paste0(prefix,"_DGF.FLR.bed"))

