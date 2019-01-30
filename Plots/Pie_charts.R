#Code to generate piecharts and barplots representing the distribution of DGF and DHS 
#author: Reyna-Llorens I, 2018.

library(RColorBrewer)
require(plyr)
require(ggplot2)
library(ggrepel)

#Set file destination
setwd("/media/iar28/7cb991c7-8087-4444-ad62-da3ae7c40dae/home/iar28/documents/DNAse-seq/Results/DGFs_wbias/Distribution_plots/DHS_dist/")

#read in files
data<-read.csv("pavis_dhs_all.csv", head=T)
data$Feature<-factor(data$Feature, levels=c("Promoter","5’UTR","Exons","Introns","3’UTR","Downstream","Intergenic"))
data$Species<-factor(data$Species, levels=c("S. italica","S. bicolor","Z. mays","B. distachyon"))

#plot
colPalette <- brewer.pal(7, "Set1")
bar <- ggplot(data, 
              aes(x = factor(1), y=Percentage,  
                  fill = data$Feature)) +
  geom_bar(colour="black",stat="identity", width=1)
bar= bar + facet_wrap(~ Species, ncol=2)
bar + coord_polar(theta = "y", direction=-1)+
  theme_void()+
  scale_fill_brewer(palette = "Set1")

#######
#################################################

Sb_WL<-read.csv("dhs_sb_wl.csv",head=T, sep="\t")
Si_WL<-read.csv("dhs_si_wl.csv",head=T, sep="\t")
Bd_WL<-read.csv("dhs_bd_wl.csv",head=T, sep="\t")
Zm_WL<-read.csv("dhs_zm_wl.csv",head=T, sep="\t")

summary(Sb_WL)

table(Sb_WL$Category[which(Sb_WL$Distance.to.TSS < 2000)])
table(Si_WL$Category[which(Sb_WL$Distance.to.TSS < 2000)])
table(Bd_WL$Category[which(Sb_WL$Distance.to.TSS < 2000)])
table(Zm_WL$Category[which(Sb_WL$Distance.to.TSS < 2000)])

f <-read.csv("1F_data_dhs.csv", head=T)
f$Species<-factor(f$Species, levels=c("S. italica","S. bicolor","Z. mays","B. distachyon"))
ggplot(f,aes(x = Species, y=value,  
             fill = Feature)) +
  geom_bar(colour="black",stat="identity",
           position="dodge",width = 0.9)+
  scale_fill_manual(values=c("gray","blue"))+
  scale_y_continuous(limits = c(0,11000), expand = c(0, 0)) +
  theme_bw()+
  theme(aspect.ratio = 1)


