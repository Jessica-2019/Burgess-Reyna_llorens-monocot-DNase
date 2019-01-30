#Script to generate circle plot for visualising the overlap between DHS, DGF and ChIPseq signals from Z. mays
#Author: Reyna-Llorens I, 2018.



setwd("/media/iar28/7cb991c7-8087-4444-ad62-da3ae7c40dae/home/iar28/documents/Zmays_ChIPseq/dnase/dnase_chipseq_AGVP3")
require(circlize)
cytoband.df = read.table("zmkar.txt", colClasses = c("character","numeric",
                                                       "numeric", "character",
                                                     "character"), sep = "\t")
head(cytoband.df)
bed = read.table("chip_seq_all_AGVPv3.bed",na.strings=NA, colClasses = c("character","numeric",
                                                "numeric", "character","numeric"), sep = "\t")

head(bed)
colnames(bed)<-c("chr","start","end","value1","value2")
bed$value1<-bed$value2

wldhs = read.table("DHS_Zm_WL.bed",na.strings=NA, colClasses = c("character","numeric",
                                                                "numeric"), sep = "\t")

head(wldhs)
wldhs$value1<-2

colnames(wldhs)<-c("chr","start","end","value1")

head(wldhs)

wldgf = read.table("zmwl_dgf.new.bed",na.strings=NA,
                   colClasses = c("character","numeric",
                   "numeric","character",
                 "numeric","character"), sep = "\t")

colnames(wldgf)<-c("chr","start","end","value1",
                   "value2","value3")

wldgfint = read.table("zmdgfchip_no_scaffolds.bed",na.strings=NA,
                   colClasses = c("character","numeric",
                                  "numeric","character",
                                  "numeric"), sep = "\t")


colnames(wldgfint)<-c("chr","start","end","value1",
                   "value2")
head(wldgfint)


circos.info()

#step 2
test2<-bed[sample(nrow(bed), 200), ]
test3<-wldhs[sample(nrow(wldhs),200),]
summary(test3)

###extend chromosomes
extend_chromosomes = function(bed, chromosome, prefix = "zoom_") {
  zoom_bed = bed[bed[[1]] %in% chromosome, , drop = FALSE]
  zoom_bed[[1]] = paste0(prefix, zoom_bed[[1]])
  rbind(bed, zoom_bed)
}
cytoband =read.cytoband(cytoband.df)
cytoband_df = cytoband$df
chromosome = cytoband$chromosome

xrange = c(cytoband$chr.len, cytoband$chr.len["1"])
xrange
normal_chr_index = 1:10
zoomed_chr_index = 11

# normalize in normal chromsomes and zoomed chromosomes separately
sector.width = c(xrange[normal_chr_index] / sum(xrange[normal_chr_index]), 
                 xrange[zoomed_chr_index] / sum(xrange[zoomed_chr_index])) 

#The extended cytoband data which is in form of a data frame is sent to circos.initializeWithIdeogram(). You can see the ideograms for chromosome 1 and 2 are zoomed (Figure 8.7).

#plot extended chromosome
circos.clear()
circos.par(start.degree = 90,
           "gap.degree" = (c(rep(c(4, 4), 5),4)))
circos.initializeWithIdeogram(extend_chromosomes(cytoband_df, "1"), plotType = c("axis", "labels"),
                              sector.width = sector.width)
circos.genomicDensity(extend_chromosomes(test3, "1"), col = c("blue"), track.height = 0.1)
circos.genomicDensity(extend_chromosomes(test2, "1"), col = c("green"), track.height = 0.1)
circos.genomicRainfall(extend_chromosomes(test2, "1"), col = c("#FF000080"),
                       pch = 16, cex = 0.4,track.height = 0.1)
circos.link("1", get.cell.meta.data("cell.xlim", sector.index = "1"),
            "zoom_1", get.cell.meta.data("cell.xlim", sector.index = "zoom_1"),
            col = "#00000020", border = NA)

circos.clear()

#plot not extended 
circos.clear()
circos.par(start.degree = 90,#"gap.degree" = (rep(c(4, 4), 5)),
           "track.height"=0.05,"track.margin"=c(0,0),"cell.padding"=c(0,0,0,0))
circos.initializeWithIdeogram(cytoband.df, chromosome.index=c(1,2,3,4,5,6,7,8,9,10),
                              plotType = c("axis", "labels"))
circos.genomicDensity(wldhs, col = c("darkblue"), track.height = 0.1)
circos.genomicDensity(wldgf, col = c("darkcyan"),
                      track.height = 0.1)
circos.genomicDensity(bed, col = c("red4"),
                       track.height = 0.1)
circos.genomicDensity(wldgfint, col = c("firebrick1"),
                      track.height = 0.1)

 circos.clear()

#only chromosome 1
circos.par("track.height" = 0.1, start.degree = 90,
           canvas.xlim = c(0, 1), canvas.ylim = c(0, 1), gap.degree = 270,
           "track.margin"=c(0,0),"cell.padding"=c(0,0,0,0))
circos.initializeWithIdeogram(cytoband.df, chromosome.index = 1,
                              plotType = c("axis","labels"))
circos.genomicDensity(wldhs, col = c("darkblue"), track.height = 0.1)
circos.genomicDensity(wldgf, col = c("darkcyan"),
                       track.height = 0.05)
circos.genomicDensity(bed, col = c("red4"),
                       track.height = 0.05)
circos.genomicDensity(wldgfint, col = c("firebrick1"),
                      track.height = 0.05)


circos.clear()
