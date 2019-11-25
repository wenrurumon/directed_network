
rm(list=ls())

#Module

rdata <- function(x){rnorm(x)}
align_mat <- function(x,y){
  #x as calculated network
  #y as reference network
  diag(x) <- 0
  diag(y) <- 0
  x <- which(as.vector(x)>0)
  y <- which(as.vector(y)>0)
  power <- mean(y%in%x)
  fdr <- 1-mean(x%in%y)
  c(POWER=power,FDR=fdr)
}
rexpr <- function(mat,nsample=500,drate=0){
  x <- matrix(0,nsample,ncol(mat))
  sel0 <- which(rowSums(mat)==0)
  for (i in sel0){
    x[,i] <- rdata(nsample)
  }
  sel <- sel0
  count.itv <- 0
  while(TRUE){
    count.itv <- count.itv+1
    sel1 <- which(rowSums(mat[,-sel0])==0)
    sel0 <- sel1[!sel1%in%sel]
    sel <- sel1
    for (i in sel0){
      x[,i] <- rowSums(x[,which(mat[i,]==1),drop=F])/2
    }
    if(count.itv>=ncol(mat)){break}
  }
  colnames(x) <- colnames(mat)
  x <- apply(x,2,minmax)
  apply(x + minmax(rnorm(length(x)))*dummy,2,minmax)
}
minmax <- function(x){
  (x-min(x))/(max(x)-min(x))
}

# Dummy data

mat <- randomDAG(30,0.1)
plot(mat)
data.mat <- unlist(mat@edgeL)
mat <- matrix(0,30,30)
for(i in 1:length(data.mat)){
  mat[data.mat[i],as.numeric(strsplit(names(data.mat[i]),'\\.')[[1]][1])] <- 1
}
dimnames(mat) <- list(paste0('g',1:ncol(mat)),paste0('g',1:ncol(mat)))

x <- apply(rexpr(mat,500),2,minmax)

#cor.matrix

cor.mat <- function(x){
  cor.mat <- apply(x,2,function(i){
    apply(x,2,function(j){
      cor.test(i,j)$p.value
    })
  })
  cor.mat <- (cor.mat <= (0.05/length(cor.mat)))
  (cor.mat+0)
}
align_mat(cor.mat(x),mat)

#TOM

library(WGCNA)
library(cluster)
options(stringsAsFactors = FALSE);


