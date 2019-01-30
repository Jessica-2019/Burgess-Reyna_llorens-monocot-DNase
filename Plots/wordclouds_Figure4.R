#Script to generate wordclouds for Figure 4
#Ivan Reyna-Llorens, 2018


require(tidyverse)
require(tm)
library(wordcloud)

setwd("/media/iar28/7cb991c7-8087-4444-ad62-da3ae7c40dae/home/iar28/documents/DNAse-seq/fimo_out/plots_fimo/")


#ARGUMENTS
args <- commandArgs(trailingOnly = TRUE)
datos <- args[1] # output from fimo Fimo_out folder *.q 
# species and celltype for tydiness 
species <-args[2] 
cellt<-args[3]

fam <- read_delim("all_annotation_noPPP2.csv", ",", col_names = T) #annotation for TF
head(fam)

length(unique(as.factor(toupper(fam$id)))) # unique TF used

colnames(fam)<-c("TF","ID","Family","source")
fam$TF<-toupper(fam$TF)
fam$ID<-toupper(fam$ID)
summary(as.factor(fam$source))
head(fam)

fam<-fam %>% unique() 
head(fam)
fam$TF<-as.factor(fam$TF)

dat<-read_delim(datos, "\t", col_names = F, na = c("", "NA"))
dat$X7 <- ifelse(grepl("MEME-1", dat$X6), dat$X4, dat$X6)
dat$X6 <- ifelse(grepl("JASPAR2018_CORE_plants_non-redundant",
                       dat$X4), dat$X5, dat$X7)

dat<-dat[,-7]
dat <-dat[,-c(4:5)] %>% filter(X1 != "NA")
colnames(dat)<- c("chr","start","end","TF")
dat$TF<-toupper(dat$TF)
dat$TF<-as.factor(dat$TF)
head(dat)
head(fam)

#frequency table counting one family member per dgf

tf_family <-data.frame(inner_join(dat, fam) %>% 
  select(chr, start,end, Family)%>%
  unique() %>% select(Family) %>%
  table())

colnames(tf_family) <-c("Var1", "Freq")
tf_family$rank <- rank(-tf_family$Freq)
tf_family<-arrange(tf_family,rank)

# Worclouds

pdf(paste0(species,"_",cellt,"_wc.pdf"))
wordcloud(rep(tf_family$Var1, tf_family$Freq),rot.per  = 0,
          scale=c(4,0.7),
          colors=brewer.pal(8,"PuOr"),random.order=FALSE)
dev.off()

#frequency table counting all family members per dgf
inner_join(dat, fam) %>% select(TF)%>%unique() %>% nrow()

tf_tf <-data.frame(inner_join(dat, fam) %>% unique() %>% select(ID) %>% table())
nrow(tf_tf)
tail(tf_tf)
colnames(tf_tf) <-c("Var1", "Freq")
tf_tf$rank <- rank(-tf_tf$Freq)
tf_tf<-arrange(tf_tf,rank)
tf_tf$species <-species
tf_tf$cell <-cellt
write.csv(tf_tf, paste0(species,"_",cellt,"_tf_ranks.known.csv"), row.names = F)



