
rm(list=ls())

source('/Users/wenrurumon/Documents/fun/rpackage/CNIF/R/sparse_2sem_final.R')
minmax <- function(x){
  (x-min(x))/(max(x)-min(x))
}
align_mat <- function(x,y){
  diag(x) <- 0
  diag(y) <- 0
  x <- which(as.vector(x)>0)
  yneg <- which(as.vector(y)==0)
  y <- which(as.vector(y)>0)
  power <- mean(y%in%x)
  fdr <- mean(x%in%yneg)
  c(POWER=power,FDR=fdr)
}
simulation <- function(seed,mp=2,Ncol1=40,Ncol2=10,Nsample=500,lambda=0.1,dummy_rate=0.01){
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
  x <- matrix(0,Nsample,ncol(mat))
  dimnames(x) <- list(paste0('s',1:Nsample),colnames(mat))
  sel0 <- which(rowSums(mat)==0)
  for (i in sel0){
    x[,i] <- minmax(runif(Nsample,0,1)) 
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
  x <- apply(x + minmax(rnorm(length(x)))*dummy_rate,2,minmax)
  #Cormat
  cor.mat <- apply(x[,1:40],2,function(i){
    apply(x[,1:40],2,function(j){
      cor.test(i,j)$p.value
    })
  })
  diag(cor.mat) <- 1
  cor.mat <- (cor.mat <= quantile(cor.mat,0.05))
  cor.mat1 <- (cor.mat+0)
  cor.mat <- apply(x,2,function(i){
    apply(x[,1:40],2,function(j){
      cor.test(i,j)$p.value
    })
  })
  diag(cor.mat) <- 1
  cor.mat <- (cor.mat <= quantile(cor.mat,0.05))
  cor.mat2 <- (cor.mat+0)
  #WGCNA
  TOM <- TOMsimilarityFromExpr(x[,1:40],power=6)
  diag(TOM) <- 0
  tom.mat1 <- (TOM>=quantile(TOM,0.95))
  TOM <- TOMsimilarityFromExpr(x,power=6)
  diag(TOM) <- 0
  tom.mat2 <- (TOM>=quantile(TOM,0.95))
  #SEM
  sem.mat1 <- sparse_2sem(x[,1:40],lambda=0.05)[[1]]
  sem.mat2 <- sparse_2sem(x,lambda=0.05)[[1]]
  #Resulting
  rbind(
    align_mat(cor.mat1,mat[1:40,1:40]+t(mat[1:40,1:40])),
    align_mat(cor.mat2,mat+t(mat)),
    align_mat(tom.mat1,mat[1:40,1:40]+t(mat[1:40,1:40])),
    align_mat(tom.mat2,mat+t(mat)),
    align_mat(sem.mat1,mat[1:40,1:40]+t(mat[1:40,1:40])),
    align_mat(sem.mat2,mat+t(mat))
  )
}
