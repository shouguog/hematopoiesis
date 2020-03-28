rm(list=ls())
#' @param corrMat correlation matrix 
matrix_zscores<-function(corrMat){
  z_row<-scale(t(corrMat))**2; 
  z_col<-scale(corrMat)**2;
  ans<-sqrt( z_row+z_col);
  cat("Zscore calculated finished....\n")
  ans;
}

#' Helper function to identify transcriptional regulators based on gene ontology annotations
#' @param species species: Hs for human or Mm for mouse

find_tfs_withGO<-function(species='human'){
  #cat("Loading gene annotations ...\n")
  library(GO.db);
    if(species=='human'){
    library(org.Hs.eg.db);
    geneSymbols<-as.list(org.Hs.egSYMBOL);
    goEntrez<-as.list(org.Hs.egGO2ALLEGS);
  }else{
    library(org.Mm.eg.db);
    geneSymbols<-as.list(org.Mm.egSYMBOL);
    goEntrez<-as.list(org.Mm.egGO2ALLEGS);
  }
  goterms<-as.list(GOTERM);
  goids<-names(goEntrez);
  onts<-lapply(goids, Ontology);
  bps<-onts[onts=='BP'];
  goids<-names(unlist(bps));
  cat("matching gene symbols and annotations")
  gobpList<-list();
  for(goid in goids){
    egs <- goEntrez[[ goid ]];
    goterm<-Term(goterms[[goid]]);
    genes<-sort(unique(as.vector(unlist( geneSymbols[egs] ))));
    gobpList[[goterm]]<-genes;
  }
  ### newHsTRs<-gobpList[['regulation of transcription, DNA-dependent']];
  regNames<-names(gobpList)[grep("regulation of transcription", names(gobpList))];
  trs<- unique(unlist(gobpList[regNames]));
  mfs<-onts[onts=='MF'];
  goidsMF<-names(unlist(mfs));
  gomfList<-list();
  for(goid in goidsMF){
    egs <- goEntrez[[ goid ]];
    goterm<-Term(goterms[[goid]]);
    genes<-sort(unique(as.vector(unlist( geneSymbols[egs] ))));
    gomfList[[goterm]]<-genes;
  }
  dbs<-gomfList[['DNA binding']];
  sort(intersect(trs, dbs));
}


#' Create global gene regulatory network
#' @param corrValue correlation
#' @param tfs list of transcriptional regulators
#' @param threshold zscore threshold to identify regulatory interactions

globalGRNFromCorr <- function(corrValue, tfs, threshold){
  tfs <- intersect(tfs, rownames(corrValue))
  corrAll <- abs(corrValue)
  zscores <- matrix_zscores(corrAll)
  zscores <- zscores[,tfs]
  grnTable <- extractRegulators(zscores, corrAll, rownames(corrValue), threshold, corrValue)
  grnTable
}
#' Helper function convert mouse gene names
#' @param genes gene names

mouseGeneNameConvert <- function(genes){
  genesCapi<-c()
  for(gene in genes){
    genesCapi<-c(genesCapi,paste0(substr(gene,1,1), tolower(substr(gene,2,100))))
  }
  genesCapi
}

#' Helper function to extract transcriptional regulators
#' @param zscores zscores
#' @param corrMatrix correlation matrix
#' @param genes genes
#' @param threshold threshold

extractRegulators <- function(zscores, corrMatrix, genes, threshold, corrValue){
#  cat(genes);cat("\t")
  #initial to speeding up
  targets <- rep('', 1e6)
  regulators <- rep('', 1e6)
  zscoresX <- rep('', 1e6)
  correlations <- rep('', 1e6)
  correlationsValue <- rep('', 1e6)
  i <- 1
  j <- 1
  for(target in genes){
    x <- zscores[target,]
    regs <- names(which(x > threshold))
    if(length(regs) > 0){
      zzs <- x[regs]
      ncount <- length(zzs)
      corrs <- corrMatrix[target, regs]
      corrsValue <- corrValue[target, regs]
      j <- i + ncount - 1
      targets[i:j] <- rep(target, ncount)
      regulators[i:j] <- regs
      zscoresX[i:j] <- zzs
      correlations[i:j] <- corrs
      correlationsValue[i:j] <- corrsValue
      i <- j + 1
    }
  }
  targets <- targets[1:j]
  regulators <- regulators[1:j]
  zscoresX <- zscoresX[1:j]
  correlations <- correlations[1:j]
  correlationsValue <- correlationsValue[1:j]
  data.frame(target=targets, reg=regulators, zscore=zscoresX, corr=correlations, correlationsValue=correlationsValue)
}

# return a iGraph object
# table of TF, TF, maybe zscores, maybe correlations
ig_tabToIgraph<-function(grnTab,  simplify=FALSE, directed=FALSE, weights=TRUE){
  tmpAns<-as.matrix(grnTab[,c("TF", "TG")]);
  regs<-as.vector(unique(grnTab[,"TF"]));
  ###cat("Length TFs:", length(regs), "\n");
  targs<-setdiff( as.vector(grnTab[,"TG"]), regs);
  ###  cat("Length TGs:", length(targs), "\n");
  myRegs<-rep("Regulator", length=length(regs));
  myTargs<-rep("Target", length=length(targs));
  types<-c(myRegs, myTargs);
  verticies<-data.frame(name=c(regs,targs), label=c(regs,targs), type=types);
  iG<-graph.data.frame(tmpAns,directed=directed,v=verticies);
  if(weights){
    #E(iG)$weight<-grnTab$weight;
    #E(iG)$weight<-as.numeric(as.character(grnTab$zscore));
    #print(range(as.numeric(as.character(grnTab$corr))))
    E(iG)$weight<-as.numeric(as.character(grnTab$corr));
    #print(range(E(iG)$weight))
  }
  if(simplify){
    #iG<-igraph::simplify(iG, remove.loops = TRUE, remove.multiple = FALSE);
    iG<-igraph::simplify(iG);
  }
  V(iG)$nEnts<-1;
  iG;
}
library('dplyr')
tfs <- find_tfs_withGO(species = "human")
library(bigSCale)

###Load the result from regNetwork
######The code from bigScale
#rm(list=ls())
#library(Seurat)
#load("../analysis/humanmouse.RData")
#human@meta.data<-human@meta.data[,c(1:3,8)]   #keep only necessary column
#mouse@meta.data<-mouse@meta.data[,c(1:3,13)] #keep only necessary column
#expr.human<-human@raw.data; rm(human)
#expr.mouse<-mouse@raw.data;remove(mouse)
##gene.names.human<-as.matrix(data.frame(V1=rownames(expr.human)))
#gene.names.mouse<-as.matrix(data.frame(V1=rownames(expr.mouse)))
#library(bigSCale)
##results.human=compute.network(expr.data = expr.human,gene.names = gene.names.human,clustering = 'recursive',speed.preset = "slow")
##save(results.human, file="networkhuman.RData")
#results.mouse=compute.network(expr.data = expr.mouse,gene.names = gene.names.mouse,clustering = 'recursive',speed.preset = "slow")
#save(results.mouse, file="networkmouse.RData")
######End The code from bigScale

load("../../regNetwork/networkhuman.RData")
grn <- globalGRNFromCorr(as.numeric(results.human$correlations), tfs, 3)
colnames(grn)[1:2]<-c("TG", "TF");
library(igraph)
ggrn<- ig_tabToIgraph(grn, simplify = TRUE,weights = TRUE)
x <- list(GRN=ggrn, GRN_table=grn, tfs=tfs)
save(x, file="ResultHumam.RData")
