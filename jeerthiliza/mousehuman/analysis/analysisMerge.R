setwd("~/work/humanmouseanalysis")
rm(list=ls())
library(Seurat)
load("humanmouse.RData")
mouse@meta.data<-mouse@meta.data[,1:3]
human@meta.data<-human@meta.data[,1:3]
speices.combined <- MergeSeurat(object1 = human, object2 = mouse, add.cell.id1 = "human", add.cell.id2 = "mouse", project = "Merge")
speices.combined
speices.combined <- FindVariableGenes(speices.combined, do.plot = F)
genes.use <- head(rownames(speices.combined@hvg.info), 1500)

speices.combined <- ScaleData(speices.combined, display.progress = F, genes.use = genes.use)
# t-SNE and Clustering
speices.combined <- RunPCA(object = speices.combined, features =genes.use)
speices.combined <- RunTSNE(speices.combined, reduction.use = "pca", dims.use = 1:20, do.fast = T)
save(speices.combined, file="speices.Merge.noalign.RData")
