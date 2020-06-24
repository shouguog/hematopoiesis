setwd("/data/gaos2/tmp/mousehuman/analysis")
rm(list=ls())
library(Seurat)
load("humanmouse.aligned.2.RData")

png("tSNE_feature_plot.png", width=2000, height=2000, res=100)
FeaturePlot(object = speices.combined, features.plot = c("GATA1", "GATA2", "CREM", "CD8A", "DPH7", "CD79A", "MLLT3", "CD37", "IL2RB"), min.cutoff = "q9", 
            cols.use = c("lightgrey", "blue"), pt.size = 0.5)
dev.off()
write.csv(speices.combined@dr$tsne@cell.embeddings, "tSNE.csv")
png("tSNE_checkData.png", width=2000, height=2000, res=100)
plot(speices.combined@dr$tsne@cell.embeddings[,1], speices.combined@dr$tsne@cell.embeddings[,2], pch=19)
dev.off()
speices.combined<-SetAllIdent(speices.combined, id = "res.0.6")
for(cluster in 14:16){
  markers <- FindConservedMarkers(speices.combined, ident.1 = cluster, grouping.var = "speice",print.bar = FALSE)
  write.csv(markers, file=paste("aligned_cluster_", cluster, ".csv",sep=""))
}

png("tSNE_celltype.png", width=2000, height=2000, res=100)
TSNEPlot(speices.combined, do.label = T, pt.size = 0.5)
dev.off()

