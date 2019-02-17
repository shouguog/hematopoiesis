setwd("H:/single_cell_RNA-seq/noncoding/PCR/Analysis/forFigures")
rm(list=ls())
library(destiny) 
library(Biobase)
library(scater)
library(scran)
#' @param new.device a logical value. If TRUE, creates a new device
#' @param bg the background color of the device
#' @param width the width of the device
rgl_init <- function(new.device = FALSE, bg = "white", width = 1640) { 
  if( new.device | rgl.cur() == 0 ) {
    rgl.open()
    par3d(windowRect = 50 + c( 0, 0, width, width ) )
    rgl.bg(color = bg )
  }
  rgl.clear(type = c("shapes", "bboxdeco"))
  rgl.viewpoint(theta = 15, phi = 20, zoom = 1)
}

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
rgl_init()
rgl.viewpoint(theta = 20, phi = 10, zoom = 1)
rgl.spheres(mapping[,1], mapping[,2], mapping[,3], r = 0.007, color = colors) 
#htmlwidgets::saveWidget(rglwidget(), "celltype.html")
writeWebGL(dir="./", filename = "celltype.html")
rgl.close()
