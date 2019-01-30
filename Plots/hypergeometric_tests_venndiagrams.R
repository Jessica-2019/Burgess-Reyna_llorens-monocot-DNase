#Code for hypergeometric test and Venn diagram generation for motif enrichment analysis between BS and WL samples
# Author: Reyna-Llorens I, 2018

require(tidyverse)
setwd("/media/iar28/7cb991c7-8087-4444-ad62-da3ae7c40dae/home/iar28/documents/DNAse-seq/fimo_out/plots_fimo/fig4g")

files = list.files(pattern="*.csv")
all<- lapply(files, read_csv) %>% bind_rows()

all$species<-c(rep("Sorghum", 870),rep("Setaria",845),rep("Maize",754))
all$cell<- c(rep("BS",417),rep("WL",453),rep("BS",432),rep("WL",413),
             rep("BS",369), rep("WL",385))
all$ID<-as.factor(toupper(all$ID))
colnames(all)<-c("tf_name","p_value","p_adj","species","cell")
fam <-read_csv("../../all_annotation_noPPP2.csv")
fam$id<-as.factor(toupper(fam$id))
fam$tf_name<-as.factor(toupper(fam$tf_name))

famall <-full_join(x=all, y=fam, by="tf_name", all=T) %>% 
  select(tf_name,p_adj,species,cell,id,Family)%>%unique()

sb <-famall %>% filter(p_adj <= 0.00099 & species =="Sorghum") 
si <-famall %>% filter(p_adj <= 0.00099 & species =="Setaria") 
zm <-famall %>% filter(p_adj <= 0.00099 & species =="Maize")
enriched<-bind_rows(sb,si,zm)
write.csv(enriched,"hypergeom_enriched.csv", row.names = F)

library(VennDiagram)

sbwl=sb%>%filter(cell=="WL")%>%select(id)%>%unique()%>%na.omit()
siwl=si%>%filter(cell=="WL")%>%select(id)%>%unique()%>%na.omit()
zmwl=zm%>%filter(cell=="WL")%>%select(id)%>%unique()%>%na.omit()

nrow(zmwl)
nrow(siwl)
nrow(sbwl)

nrow(zmbs)
nrow(sibs)
nrow(sbbs)

rbind(sbwl,siwl,zmwl) %>%unique() %>% nrow #133 motifs enriched in wl
rbind(sbbs,sibs,zmbs) %>%unique() %>% nrow #106 in bs
sbbs=sb%>%filter(cell=="BS")%>%select(id)%>%unique()%>%na.omit()
sibs=si%>%filter(cell=="BS")%>%select(id)%>%unique()%>%na.omit()
zmbs=zm%>%filter(cell=="BS")%>%select(id)%>%unique()%>%na.omit()

library(VennDiagram)
pdf("venn_BS.pdf")
venn.plot <- venn.diagram(list(sibs$id, sbbs$id, zmbs$id),
                          NULL,
                          fill=c("orange","orange","orange"),
                          lwd=c(1,1,1),
                          alpha=c(0.3,0.2,0.1),
                          cex = 2,
                          cat.fontface=4,
                          category.names = c("Setaria" , "Sorghum" , "Maize"),
                          cat.fontfamily = "sans",
                          rotation = 1
)
grid.draw(venn.plot)
dev.off()

pdf("venn_WL.pdf")
venn.plot <- venn.diagram(list(siwl$id, sbwl$id, zmwl$id),
                          NULL,
                          fill=c("blue","blue","blue"),
                          lwd=c(1,1,1),
                          alpha=c(0.3,0.2,0.1),
                          cex = 2,
                          cat.fontface=4,
                          category.names = c("Setaria" , "Sorghum" , "Maize"),
                          cat.fontfamily = "sans",
                          rotation = 1
)
grid.draw(venn.plot)
dev.off()



b<-read_csv("hypergeom_enriched.csv")
c<-read_csv("../tf_family.csv", col_names = F)
colnames(c)<-c("tf_name","fam","gene")
head(c)
c$tf_name<-toupper(c$tf_name)

bc <-left_join(x=b, y=c, by="tf_name", all=T)%>% unique() 
tail(bc)
write.csv(bc,"hyp.tmp2.csv",row.names = F)

head(bc)

