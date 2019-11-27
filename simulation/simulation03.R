
rm(list=ls())

library(pcalg)
library(WGCNA)
library(dplyr)
library(RcppArmadillo)
setwd('e://huzixin//')
source('sparse_2sem.R')
source('CNIF.R')
source('local_cnif_macro.R')
Rcpp::sourceCpp('initial_sem.cpp')
Rcpp::sourceCpp('score_function_regression.cpp')
Rcpp::sourceCpp('simple_cycle.cpp')

minmax <- function(x){
  (x-min(x))/(max(x)-min(x))
}
align_mat <- function(x,y,rebase=T){
  diag(x) <- 0
  diag(y) <- 0
  if(rebase){
    y <- y+t(y)
    x <- x+t(x)
  }
  xneg <- which(as.vector(x)==0)
  xpos <- which(as.vector(x)>0)
  yneg <- which(as.vector(y)==0)
  ypos <- which(as.vector(y)>0)
  power <- mean(ypos%in%xpos)
  fdr <- mean(yneg%in%xpos)
  c(
    # x=length(x),y=length(y),
    # lap=sum(x%in%y),all=length(x)+length(xneg),
    POWER=power,FDR=fdr)
}
align <- function(x,y){
  # x <- cor.mat2
  # y <- mat
  # xneg <- which(as.vector(x)==0)
  # xneg2 <- which(as.vector(x+t(x))==0)
  # xpos <- which(as.vector(x)>0)
  # xpos2 <- which(as.vector(x+t(x))>0)
  # yneg <- which(as.vector(y)==0)
  # yneg2 <- which(as.vector(y+t(y))==0)
  # ypos <- which(as.vector(y)>0)
  # ypos2 <- which(as.vector(y+t(y))>0)
  # mean(ypos%in%xpos)
  # mean(ypos2%in%xpos2)
  # power1 <- sum((x>0)*(y>0))/sum(y>0)
  # power2 <- sum(((x+t(x))>0)*((y+t(y))>0))/sum((y+t(y))>0)
  rlt1 <- align_mat(x,y,T)
  rlt2 <- align_mat(x,y,F)
  c(rlt1,rlt2)
}
simu <- function(seed,mp=2,Ncol1=40,Ncol2=10,Nsample=500,lambda=0.1,dummy_rate=0.01){
  #setup
  # seed <- 1000
  # mp <- 2
  # Ncol1 <- 40
  # Ncol2 <- 10
  # Nsample <- 500
  # lambda <- 0.1
  # dummy_rate <- 0.1
  #Module
  set.seed(seed)
  #Random DAG
  mat <- randomDAG(Ncol1,lambda)
  data.mat <- unlist(mat@edgeL)
  mat <- matrix(0,Ncol1,Ncol1)
  for(i in 1:length(data.mat)){
    mat[data.mat[i],as.numeric(strsplit(names(data.mat[i]),'\\.')[[1]][1])] <- 1
  }
  for(i in 1:nrow(mat)){
    if(sum(mat[i,])>mp){
      mat[i,sample(which(mat[i,]>0),sum(mat[i,])-mp)] <- 0
    }
  }
  mat11 <- mat
  dimnames(mat11) <- list(paste0('e',1:Ncol1),paste0('e',1:Ncol1))
  mat12 <- matrix(0,Ncol1,Ncol2)
  dimnames(mat12) <- list(paste0('e',1:Ncol1),paste0('g',1:Ncol2))
  mat12[rowSums(mat11)==0,] <- apply(mat12[rowSums(mat11)==0,],1,function(x){
    x[sample(1:length(x),sample(c(1,1,1,2,2,3))[1])] <- 1
    x
  })
  mat12[rowSums(mat11)>0,] <- apply(mat12[rowSums(mat11)>0,],1,function(x){
    x[sample(1:length(x),sample(c(0,0,0,0,1,1,1,2,3))[1])] <- 1
    x
  })
  mat12[which(rowSums(cbind(mat11,mat12))==0),1] <- 1
  mat <- cbind(mat11,mat12)
  mat <- rbind(mat,matrix(0,Ncol2,Ncol1+Ncol2))
  rownames(mat) <- colnames(mat)
  #Random Data
  x <- matrix(rnorm(Nsample*ncol(mat),0,dummy_rate),Nsample,ncol(mat))
  dimnames(x) <- list(paste0('s',1:Nsample),colnames(mat))
  sel0 <- which(rowSums(mat)==0)
  for (i in sel0){
    # x[,i] <- minmax(runif(Nsample,0,1)) 
    x[,i] <- x[,i]+sample(c(0,1,2),size=Nsample,replace=T) %>% scale
  }
  sel <- sel0
  while(TRUE){
    sel1 <- which(rowSums(mat[,-sel,drop=F])==0)#nodes wo any mp after sel
    sel0 <- sel1[!sel1%in%sel] #nodes to be dummied
    sel <- sel1 #nodes dummied
    if(length(sel0)==0){break}
    for (i in sel0){
      x1 <- x[,which(mat[i,]==1),drop=F]
      c1 <- matrix(runif(sum(mat[i,]==1),1,2),nrow=sum(mat[i,]==1))
      x[,i] <- scale(x1 %*% c1) + x[,i]
    }
  }
  colnames(x) <- colnames(mat)
  #Cormat
  cor.mat <- apply(x[,1:40],2,function(i){
    apply(x[,1:40],2,function(j){
      cor.test(i,j)$p.value
    })
  })
  diag(cor.mat) <- 1
  cor.mat <- (cor.mat <= quantile(cor.mat,0.1))
  cor.mat1 <- (cor.mat+0)
  cor.mat <- apply(x,2,function(i){
    apply(x,2,function(j){
      cor.test(i,j)$p.value
    })
  })
  diag(cor.mat) <- 1
  cor.mat <- (cor.mat <= quantile(cor.mat,0.1))
  cor.mat2 <- (cor.mat+0)
  cor.mat2[1:40,1:40] <- cor.mat2[1:40,1:40] + cor.mat1
  rlt1 <- list(
    cor.mat1 = align(cor.mat1,mat[1:40,1:40]),
    cor.mat2 = align(cor.mat2,mat)
  )
  #WGCNA
  TOM <- TOMsimilarityFromExpr(x[,1:40],power=6)
  diag(TOM) <- 0
  tom.mat1 <- (TOM>=quantile(TOM,0.9))
  TOM <- TOMsimilarityFromExpr(x,power=6)
  diag(TOM) <- 0
  tom.mat2 <- (TOM>=quantile(TOM,0.9))
  tom.mat2[-1:-Ncol1] <- 0
  tom.mat2[1:40,1:40] <- tom.mat2[1:40,1:40]+tom.mat1
  rlt2 <- list(
    tom.mat1 = align(tom.mat1,mat[1:40,1:40]),
    tom.mat2 = align(tom.mat2,mat)
  )
  #SEM
  sem.mat1 <- sparse_2sem(x[,1:40],lambda=0.1)[[1]]
  sem.mat2 <- sparse_2sem(Y=x[,1:40],Y.fixed=sem.mat1,
                          X=x[,-1:-40],lambda=0.001)[[1]]
  sem.mat2 <- rbind(sem.mat2,matrix(0,Ncol2,Ncol1+Ncol2))
  rlt3 <- list(
    sem.mat1 = align(sem.mat1,mat[1:40,1:40]),
    sem.mat2 = align(sem.mat2,mat)
  )
  sem.mat3 <- CNIF(x[,1:40],init.adj=sem.mat1,max_parent=mp)
  sem.mat4 <- sparse_2sem(Y=x[,1:40],Y.fixed=sem.mat3,
                          X=x[,-1:-40],lambda=0)[[1]]
  sem.mat4[,1:40] <- sem.mat1-sem.mat3
  sem.mat4 <- rbind(sem.mat4,matrix(0,Ncol2,Ncol1+Ncol2))
  rlt4 <- list(
    sem.mat3 = align(sem.mat1-sem.mat3,mat[1:40,1:40]),
    sem.mat4 = align(sem.mat4,mat)
  )
  #Resulting
  do.call(rbind,do.call(c,list(rlt1,rlt2,rlt3,rlt4)))
}

test100 <- lapply(1:100,simu,Nsample=100) %>% dst
test300 <- lapply(1:100,simu,Nsample=300) %>% dst
test500 <- lapply(1:100,simu,Nsample=500) %>% dst
test1000 <- lapply(1:100,simu,Nsample=1000) %>% dst

dst <- function(x){
  p1 <- sapply(x,function(x){x[,1]})
  p2 <- sapply(x,function(x){x[,3]})
  f1 <- sapply(x,function(x){x[,2]})
  f2 <- sapply(x,function(x){x[,4]})
  rlt <- cbind(
    poewr1 = apply(p1,1,mean),
    FDR1 = apply(f1,1,mean),
    power2 = apply(p2,1,mean),
    FDR1 = apply(f2,1,mean)
  )
  rownames(rlt) <- c(
    'correlation','correlation_wgeno',
    'coexpression','coexpression_wgeno',
    'SEM','SEM_wgeno',
    'CNIF','CNIF_wgeno'
  )
  list(rlt=x,summary=rlt)
}
