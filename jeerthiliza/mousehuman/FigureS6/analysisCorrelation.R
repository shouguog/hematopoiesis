setwd("/data/gaos2/tmp/mousehuman/FigureS6")
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
dist<-1-cor(t(data.combined))

clusters <- hclust(as.dist(dist))
png("hcl_mouse_human_cluster.png", width=1000, height=1000, res=150)
plot(clusters, main = "conserved clustering", xlab = "clusters")
dev.off()

