#setwd("/data/gaos2/tmp/mousehuman/FigureS6")
rm(list=ls())
library(Seurat)
load("../analysis/humanmouse.RData")
human <- NormalizeData(human)
mouse <- NormalizeData(mouse)

# Gene selection for input to CCA
human <- FindVariableGenes(human, do.plot = F)
mouse <- FindVariableGenes(mouse, do.plot = F)
g.1 <- head(rownames(human@hvg.info), 4000)
g.2 <- head(rownames(mouse@hvg.info), 4000)
genes.use <- unique(c(g.1, g.2))
human <- ScaleData(human, display.progress = F, genes.use = genes.use)
mouse <- ScaleData(mouse, display.progress = F, genes.use = genes.use)

mouseData<-as.data.frame(t(mouse@scale.data))
mouseData$celltype<-mouse@meta.data$celltype
mouseClusterData<-aggregate(mouseData, list(mouseData$celltype),mean)
rownames(mouseClusterData)<-paste("mouse", mouseClusterData$Group.1, sep="_")
mouseClusterData<-mouseClusterData[, c(-1, -5832)]
mouseClusterData<-mouseClusterData[rownames(mouseClusterData)!="mouse_NotypeDefination",]

humanData<-as.data.frame(t(human@scale.data))
humanData$celltype<-human@meta.data$celltype
humanClusterData<-aggregate(humanData, list(humanData$celltype),mean)
rownames(humanClusterData)<-paste("human", humanClusterData$Group.1, sep="_")
humanClusterData<-humanClusterData[, c(-1, -5832)]

data.combined<-rbind(humanClusterData, mouseClusterData)

source("../referencenew/specTau/specTau.R")
load("../referencenew/specTau/qnExp.rda")
#head(qnExp)
tauExp <- calcTau(qnExp)
#head(tauExp)
###
library(preprocessCore)
###scale to 0-1
data.combined<-(data.combined-min(data.combined))/(max(data.combined)-min(data.combined))
data_norm <- as.data.frame(normalize.quantiles(as.matrix(t(data.combined)), copy = TRUE))
rownames(data_norm)<-colnames(data.combined)
colnames(data_norm)<-rownames(data.combined)
rm(g.1,g.2,genes.use,human,humanClusterData,humanData,mouse,mouseClusterData,mouseData)
tauValue <- calcTau(data_norm)
write.csv(tauValue, "ourDataTau.csv", quote = FALSE)
