rm(list=ls())
library(SingleCellExperiment)
library(scmap)
rm(list=ls())
maxNum<-10000
load("../database/human.RObj")
load("../database/mouseDirty_batch1.RObj")
rownames(mouseDirty@raw.data)<-toupper(rownames(mouseDirty@raw.data))
geneovelap<-intersect(rownames(GATA2_HD@raw.data), rownames(mouseDirty@raw.data))
rm(mouseDirty)  #save space
counts<-as.matrix(GATA2_HD@raw.data[geneovelap,1:maxNum])
celltype1<-GATA2_HD@meta.data$celltype[1:maxNum]
rm(GATA2_HD)
scaleRow<-colSums(counts)
countsNorm<-counts
for(ii in 1:length(scaleRow)){
  countsNorm[,ii]<-counts[,ii]/scaleRow[ii]*1000000
}
countsNormLog<-log10(countsNorm+0.2)

h <- SingleCellExperiment(assays=list(counts=counts, normcounts=countsNorm, logcounts=countsNormLog))
rowData(h)$feature_symbol<-geneovelap
colData(h)$cell_type1<-celltype1

load("../database/mouseDirty_batch1.RObj")
rownames(mouseDirty@raw.data)<-toupper(rownames(mouseDirty@raw.data))
counts<-as.matrix(mouseDirty@raw.data[geneovelap,1:maxNum])
celltype1<-mouseDirty@meta.data$celltype[1:maxNum]
counts<-counts[, celltype1!="NotypeDefination"]
celltype1<-celltype1[celltype1!="NotypeDefination"]
rm(mouseDirty)
scaleRow<-colSums(counts)
countsNorm<-counts
for(ii in 1:length(scaleRow)){
  countsNorm[,ii]<-counts[,ii]/scaleRow[ii]*1000000
}
countsNormLog<-log10(countsNorm+0.2)
m <- SingleCellExperiment(assays=list(counts=counts, normcounts=countsNorm, logcounts=countsNormLog))
rowData(m)$feature_symbol<-geneovelap
colData(m)$cell_type1<-celltype1
rm(counts, countsNorm, countsNormLog, scaleRow, celltype1)

####Use the mouse as reference
m <- selectFeatures(m)
m <- indexCluster(m)
scmapCluster_results <- scmapCluster(projection = h,index_list = list(metadata(m)$scmap_cluster_index),threshold = 0.1)
#Percentage of assigned cells:
length(which(scmapCluster_results$scmap_cluster_labs != "unassigned")) / ncol(h) * 100
#library(irr)
#kappa2(data.frame(scmapCluster_results$scmap_cluster_labs, h$cell_type1)[scmapCluster_results$scmap_cluster_labs,])
plot(getSankey(colData(h)$cell_type1,scmapCluster_results$scmap_cluster_labs,plot_height = 800, plot_width = 800))

####Use the human as reference
h <- selectFeatures(h)
h <- indexCluster(h)
scmapCluster_results <- scmapCluster(projection = m,index_list = list(metadata(h)$scmap_cluster_index),threshold = 0.1)
#Percentage of assigned cells:
length(which(scmapCluster_results$scmap_cluster_labs != "unassigned")) / ncol(m) * 100
#library(irr)
#kappa2(data.frame(scmapCluster_results$scmap_cluster_labs, h$cell_type1)[scmapCluster_results$scmap_cluster_labs,])
plot(getSankey(colData(m)$cell_type1,scmapCluster_results$scmap_cluster_labs
               ,plot_height = 800, plot_width = 800, colors=c("blue", "red", "brown", "orange", "green", "violet")))

plot(getSankey(colData(m)$cell_type1,scmapCluster_results$scmap_cluster_labs,plot_height = 800, plot_width = 800))
