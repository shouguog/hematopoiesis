setwd("/data/gaos2/tmp/mousehuman/analysis")
rm(list=ls())
library(Seurat)
load("humanmouse.RData")
human@meta.data<-human@meta.data[,c(1:3,8)]   #keep only necessary column
mouse@meta.data<-mouse@meta.data[,c(1:3,13)] #keep only necessary column
# Set up human object
human@meta.data$speice <- "human"
human <- NormalizeData(human)
# Set up mouse object
mouse@meta.data$speice <- "mouse"
mouse <- NormalizeData(mouse)

# Gene selection for input to CCA
human <- FindVariableGenes(human, do.plot = F)
mouse <- FindVariableGenes(mouse, do.plot = F)
g.1 <- head(rownames(human@hvg.info), 1000)
g.2 <- head(rownames(mouse@hvg.info), 1000)
genes.use <- unique(c(g.1, g.2))
human <- ScaleData(human, display.progress = F, genes.use = genes.use)
mouse <- ScaleData(mouse, display.progress = F, genes.use = genes.use)

speices.combined <- RunCCA(human, mouse, genes.use = genes.use, num.cc = 30)
# visualize results of CCA plot CC1 versus CC2 and look at a violin plot
p1 <- DimPlot(object = speices.combined, reduction.use = "cca", group.by = "speice", pt.size = 0.5, do.return = TRUE)
p2 <- VlnPlot(object = speices.combined, features.plot = "CC1", group.by = "speice", do.return = TRUE)
plot_grid(p1, p2)
speices.combined <- AlignSubspace(speices.combined, reduction.type = "cca", grouping.var = "speice",dims.align = 1:20)
p1 <- VlnPlot(object = speices.combined, features.plot = "ACC1", group.by = "speice", do.return = TRUE)
p2 <- VlnPlot(object = speices.combined, features.plot = "ACC2", group.by = "speice", do.return = TRUE)
plot_grid(p1, p2)
rm(human.data, mouse.data)
save(list=ls(), file="humanmouse.aligned.RData")
# t-SNE and Clustering
speices.combined <- RunTSNE(speices.combined, reduction.use = "cca.aligned", dims.use = 1:20, do.fast = T)
speices.combined <- FindClusters(speices.combined, reduction.type = "cca.aligned", resolution = 0.6, dims.use = 1:20)
# Visualization
p1 <- TSNEPlot(speices.combined, do.return = T, pt.size = 0.5, group.by = "speice")
p2 <- TSNEPlot(speices.combined, do.label = T, do.return = T, pt.size = 0.5)
png("tSNE_group_mouse_human.png", width=2000, height=1000, res=100)
plot_grid(p1, p2)
dev.off()
save(list=ls(), file="humanmouse.aligned.2.RData")

