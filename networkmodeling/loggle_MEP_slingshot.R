##R4.1.0
#setwd("/data/gaos2/tmp/2021-11/networkBiology/analysis_HDmito10/loggle")
rm(list=ls())
load("pd_scaleData.RData")
load("sling_monocle2.RData")
pd$Pseudotime_sling<-pdplusD3[rownames(pd), "PseudotimeSlingShot"]
library(loggle)
pdInterest<-pd[pd$celltype %in% c("HSC", "MEP"),]
pdInterest<-pdInterest[order(pdInterest$Pseudotime_sling),]
X <- scaleData[, rownames(pdInterest)]
###Let us bin
Y<-X
binNum<- floor(dim(Y)[2]/30)
X<-matrix(1, nrow = dim(Y)[1], ncol =binNum); rownames(X)<-rownames(Y)
for(ii in 1:binNum){
  X[,ii]<-rowMeans(Y[, -29:0+ii*30])
}
###End bin
X<-X[names(sort(apply(X, 1, sd), decreasing = TRUE))[1:80],]###Some genes have sd of 0 error
p <- nrow(X)  ## number of variables
## positions of time points to estimate graphs
pos <- round(seq(0.05, 0.85, length=25)*(ncol(X)-1)+1)
K <- length(pos)
## estimate time-varying graphs and conduct model selection via cross-validation
## num.thread can be set as large as number of cores on a multi-core machine 
## (however when p is large, memory overflow should also be taken caution of)
ts <- proc.time()
result <- loggle.cv(X, pos, h.list = c(0.15, 0.2, 0.25, 0.3), d.list = c(0, 0.01, 0.05, 0.1, 0.2, 0.3, 1), detrend = FALSE, 
                    lambda.list = c(0.15, 0.2, 0.25, 0.3), fit.type = "pseudo", cv.vote.thres = 1, num.thread = 1, cv.fold = 4)
te <- proc.time()
## optimal values of h, d and lambda, and number of selected edges at each time point
select.result <- result$cv.select.result
###Let us check the result
sumData<-cbind("time" = seq(0.05, 0.85, length=25), "h.opt" = rep(select.result$h.opt, K), "d.opt" = select.result$d.opt,
               "lambda.opt" = select.result$lambda.opt, "edge.num.opt" = select.result$edge.num.opt)
write.csv(sumData, file = "sumData_HSCtoMEP_rep1_slingshot.csv", quote = FALSE, row.names = FALSE)
save(result, X, file = "result_HSCtoMEP_rep1_slingshot.RData")
