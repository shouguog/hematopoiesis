setwd("~/work/humanmouseanalysis")
rm(list=ls())
library(Seurat)
load("speices.Merge.noalign.RData")
matadata<-speices.combined@meta.data
colors<-rep("blue", dim(matadata)[1])
colors[grepl("human_", rownames(matadata))]<-"red"
png("beforealigntSNE.png", width = 2000, height=2000, res = 100)
plot(speices.combined@dr$tsne@cell.embeddings[,1], speices.combined@dr$tsne@cell.embeddings[,2], col=colors, pch=16)
dev.off()
