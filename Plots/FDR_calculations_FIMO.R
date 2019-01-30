#Calculate FDR after running FIMO using previously defined motifs from Z.mays

require(tidyverse)

setwd("~/Desktop")

a <- read_tsv("fimo_zmwl.tsv")
a$fdr <-p.adjust(a$`p-value`, method = "BH", n = length(a$`p-value`))
a$motif_id<-as.factor(a$motif_id)
a<-a%>%select(motif_id, `p-value`,`q-value`, matched_sequence,fdr)%>%na.omit()%>% filter(`fdr` <= 0.001) 
summary(as.factor(a$motif_id))

