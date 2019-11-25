
rm(list=ls())
setwd('E:\\sample')

source("sparse_2sem_final.r")
source('local_cnif_macro.R')
source('CNIF.R')
sourceCpp("score_function_regression.cpp")
sourceCpp("simple_cycle.cpp")
sourceCpp("initial_sem.cpp")
library(pcalg)
library(WGCNA)

rdata <- function(x){rnorm(x)}
align_mat <- function(x,y){
  diag(x) <- 0
  diag(y) <- 0
  x <- which(as.vector(x)>0)
  y <- which(as.vector(y)>0)
  power <- mean(y%in%x)
  fdr <- 1-mean(x%in%y)
  c(POWER=power,FDR=fdr)
}
rexpr <- function(mat,nsample=500,dummy=0){
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

simulation <- function(seed,mp=2,Ncol=50,Nsample=500,lambda=0.1,dummy_rate=0.1){
  set.seed(seed)
  mat <- randomDAG(Ncol,lambda)
  data.mat <- unlist(mat@edgeL)
  mat <- matrix(0,Ncol,Ncol)
  for(i in 1:length(data.mat)){
    mat[data.mat[i],as.numeric(strsplit(names(data.mat[i]),'\\.')[[1]][1])] <- 1
  }
  for(i in 1:nrow(mat)){
    if(sum(mat[i,])>mp){
      mat[i,sample(which(mat[i,]>0),sum(mat[i,])-mp)] <- 0
    }
  }
  dimnames(mat) <- list(paste0('g',1:ncol(mat)),paste0('g',1:ncol(mat)))
  x <- apply(rexpr(mat,Nsample,dummy_rate),2,minmax)
  #cor.mat
  cor.mat <- apply(x,2,function(i){
    apply(x,2,function(j){
      cor.test(i,j)$p.value
    })
  })
  diag(cor.mat) <- 1
  cor.mat <- (cor.mat <= quantile(cor.mat,0.05))
  cor.mat <- (cor.mat+0)
  #sem.mat
  sem.mat <- sparse_2sem(x,lambda=0.1)[[1]] 
  sem.cnif <- CNIF(data = x,init.adj = sem.mat,max_parent = mp)
  #tom
  TOM <- TOMsimilarityFromExpr(x,power=6)
  dimnames(TOM) <- list(colnames(x),colnames(x))
  diag(TOM) <- 0
  cenet <- (TOM>=quantile(TOM,0.95));sum(cenet)
  #result
  rlt <- list(
    cor = align_mat(cor.mat,(mat+t(mat))),
    sem = align_mat(sem.mat,mat+t(mat)),
    cnif = align_mat(sem.cnif+t(sem.cnif),mat+t(mat)),
    tom = align_mat(cenet,(mat+t(mat))),
    cnif2 = align_mat(sem.cnif,mat)
  )
  rlt
}

s_100 <- lapply(1:1000,function(i){
  print(paste('Simulation',i))
  simulation(i,2,50,100,0.1,0.1)
})
s_300 <- lapply(1:1000,function(i){
  print(paste('Simulation',i))
  simulation(i,2,50,300,0.1,0.1)
})
s_500 <- lapply(1:1000,function(i){
  print(paste('Simulation',i))
  simulation(i,2,50,500,0.1,0.1)
})
s_1000 <- lapply(1:1000,function(i){
  print(paste('Simulation',i))
  simulation(i,2,50,1000,0.1,0.1)
})

