setwd("H:/single_cell_RNA-seq/noncoding/PCR/Analysis/forFigures")
rm(list=ls())
library(destiny) 
library(Biobase)
library(scater)
library(scran)
load("../../../../check_delong/allpatients/count/delongcelland391cell/destinty/difusionMapResult.RData")
library(rgl)
# Make a scatter plot
colors<-rep("grey",391)
for(ii in 1:391){
  #if(IDs[ii, "NP38"]=="38N"){colors[ii]="blue"}
  if(IDs[ii, "celltype"]=="MEP"){colors[ii]="violetred1"}
  if(IDs[ii, "celltype"]=="ETP"){colors[ii]="darkgoldenrod4"}
  if(IDs[ii, "celltype"]=="yHSC"){colors[ii]="brown3"}
  if(IDs[ii, "celltype"]=="MLP"){colors[ii]="darkgreen"}
  if(IDs[ii, "celltype"]=="GMP"){colors[ii]="cyan3"}
  if(IDs[ii, "celltype"]=="ProB"){colors[ii]="dodgerblue1"}
}
IDsname<-gsub("_subread.sorted", "", gsub("Sample_", "", gsub(".bam", "", rownames(IDs))))
IDsname<-gsub("H3.0", "H3.", gsub("H4.0", "H4.", gsub("H5_", "H5.", IDsname)))
rownames(IDs)<-IDsname
plot3d(mapping[,1], mapping[,2], mapping[,3], size = 2, col = colors,
       type="s",xlab="DM1",ylab = "DM2",zlab="DM3")
pp<-dget("PCR.R")
par3d(pp)
###manually rotate to get good view
#pp<-par3d(no.readonly = TRUE)
#dput(pp, file="PCR.R", control = "all")

writeWebGL(dir="./", filename = "celltype.html", width=600, height = 600)
####Get PCR data
rownames(mapping)<-IDsname
dataComNormByGene<-read.csv("../dataComNormByGeneWithCelltypeWithGeneAnno_nobulk.csv", header=T)
dataComNormByGeneH3<-dataComNormByGene[grepl("H3", dataComNormByGene$X),]###Final
dataComNormByGeneH3<-dataComNormByGene[grepl("H3", dataComNormByGene$X),]
dataComNormByGeneH4<-dataComNormByGene[grepl("H4", dataComNormByGene$X),]
dataComNormByGeneH5<-dataComNormByGene[grepl("H5", dataComNormByGene$X),]
dataComNormByGeneZhao<-dataComNormByGene[grepl("Zhao", dataComNormByGene$X),]
dataComNormByGeneHD<-rbind(dataComNormByGeneH3, dataComNormByGeneH4, dataComNormByGeneH5, dataComNormByGeneZhao)
###Get Overlap
dataComNormByGeneHD<-dataComNormByGeneHD[dataComNormByGeneHD$X %in%IDsname,]
###For gene filtering
genes<-read.csv("geneListUsed.csv", row.names = 1)
#genesMEP<-rownames(genes[genes$type=="MEP" & genes$use=="yes",])
genesMEP<-rownames(genes[genes$type=="MEP",])
genes38N<-rownames(genes[genes$type=="38N",])
genesLMP<-rownames(genes[genes$type=="LMP",])
dataComNormByBranch<-as.data.frame(matrix(1, nrow=dim(dataComNormByGeneHD)[1], ncol=3))
colnames(dataComNormByBranch)<-c("MEP", "LMP", "38N")
rownames(dataComNormByBranch)<-dataComNormByGeneHD$X
dataComNormByBranch$MEP<-as.numeric(rowMeans(dataComNormByGeneHD[,genesMEP]))
dataComNormByBranch[, "38N"]<-as.numeric(rowMeans(dataComNormByGeneHD[,genes38N]))
dataComNormByBranch$LMP<-as.numeric(rowMeans(dataComNormByGeneHD[,genesLMP]))
mapping<-mapping[rownames(dataComNormByBranch),]
###Let us try GATA1
colors<-c()
for(ii in 1:dim(mapping)[1]){
  if(dataComNormByGeneHD[ii, "GATA1_annoted_MEP"]<6){
    colors<-c(colors,"red")
  }else{
    colors<-c(colors,"black")
  }
}
plot3d(mapping[,1], mapping[,2], mapping[,3], size = 2, col = colors,
       type="s",xlab="DM1",ylab = "DM2",zlab="DM3")
writeWebGL(dir="./", filename = "PCRGATA1.html", width=600, height = 600)

###Let us try MEIS1
colors<-c()
for(ii in 1:dim(mapping)[1]){
  if(dataComNormByGeneHD[ii, "MEIS1_annoted_yHSC"]<6){
    colors<-c(colors,"red")
  }else{
    colors<-c(colors,"black")
  }
}
plot3d(mapping[,1], mapping[,2], mapping[,3], size = 2, col = colors,
       type="s",xlab="DM1",ylab = "DM2",zlab="DM3")
writeWebGL(dir="./", filename = "PCR_MEIS.html", width=600, height = 600)

###Let us try CD79A_annoted_GMP
colors<-c()
for(ii in 1:dim(mapping)[1]){
  if(dataComNormByGeneHD[ii, "CD79A_annoted_GMP"]<6){
    colors<-c(colors,"red")
  }else{
    colors<-c(colors,"black")
  }
}
plot3d(mapping[,1], mapping[,2], mapping[,3], size = 2, col = colors,
       type="s",xlab="DM1",ylab = "DM2",zlab="DM3")
writeWebGL(dir="./", filename = "PCR_CD79A.html", width=600, height = 600)

###Let us try CXCR4_annoted_ETP
colors<-c()
for(ii in 1:dim(mapping)[1]){
  if(dataComNormByGeneHD[ii, "CXCR4_annoted_ETP"]<6){
    colors<-c(colors,"red")
  }else{
    colors<-c(colors,"black")
  }
}
plot3d(mapping[,1], mapping[,2], mapping[,3], size = 2, col = colors,
       type="s",xlab="DM1",ylab = "DM2",zlab="DM3")
writeWebGL(dir="./", filename = "PCR_CXCR4.html", width=600, height = 600)


