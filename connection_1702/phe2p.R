
rm(list=ls())
setwd('C:\\Users\\zhu2\\Documents\\getpathway')
load("~/getpathway/gene39761.rda.RData")
load("C:/Users/zhu2/Documents/getpathway/ptwmap2.rda")

qpca <- function(A,rank=0,ifscale=TRUE){
  if(ifscale){A <- scale(as.matrix(A))[,]}
  A.svd <- svd(A)
  if(rank==0){
    d <- A.svd$d
  } else {
    d <- A.svd$d-A.svd$d[min(rank+1,nrow(A),ncol(A))]
  }
  d <- d[d > 1e-8]
  r <- length(d)
  prop <- d^2; info <- sum(prop)/sum(A.svd$d^2);prop <- cumsum(prop/sum(prop))
  d <- diag(d,length(d),length(d))
  u <- A.svd$u[,1:r,drop=F]
  v <- A.svd$v[,1:r,drop=F]
  x <- u%*%sqrt(d)
  y <- sqrt(d)%*%t(v)
  z <- x %*% y
  rlt <- list(rank=r,X=x,Y=y,Z=x%*%y,prop=prop,info=info)
  return(rlt)
}
pca <- function(X){
  X <- scale(as.matrix(X))
  m = nrow(X)
  n = ncol(X)
  Xeigen <- svd(X)
  value <- (Xeigen$d)^2/m
  value <- cumsum(value/sum(value))
  score <- X %*% Xeigen$v
  mat <- Xeigen$v
  list(score=score,prop=value,mat=mat)
}

expinpath <- expinpath[-3]
inputs <- lapply(unique(pathlist[,1]),function(grpi){
  paths <- pathlist[pathlist[,1]==grpi,2]
  input <- expinpath[names(expinpath)%in%paths]
  print(list(groupname=grpi,paths=names(input)))
  input
})
names(inputs) <- unique(pathlist[,1])
inputs <- inputs[sapply(inputs,length)>0]

pcas <- lapply(inputs,function(x){
  x <- do.call(cbind,x)
  x <- qpca(x,which(pca(x)$prop>=.9)[1])
  x$X[,1:which(x$prop>=0.9)[1],drop=F]
})

################################
# Load 2
################################

rm(list=ls()[-grep('pcas',ls())])
library(igraph)
library(WGCNA)
library(grplasso)
library(slam)

setwd('C:\\Users\\zhu2\\Documents\\signaling\\')
source('codes/sparse_2sem_final.R')
source('codes/local_cnif_macro.R')
source('codes/flm_and_cca.R')
source('codes/CNIF.R')
sourceCpp("codes/score_function_regression.cpp")
sourceCpp("codes/simple_cycle.cpp")
sourceCpp("codes/initial_sem.cpp")
source('codes/CNIF_grouplasso2.R')
setwd('C:\\Users\\zhu2\\Documents\\getpathway\\model20170215')
load('phenet.rda')
load('rlt_p2pcp.rda')
load('rlt_p2pinp.rda')

plotnet <- function(x,mode='undirected'){
  diag(x) <- 0
  plot(graph_from_adjacency_matrix(t(x),mode=mode),
       edge.arrow.size=.1,
       vertex.size=3,
       vertex.label.cex=1,
       edge.width=.1)
}
plotnet2 <- function(phe_pg,kick=FALSE,thres=0.8,main=""){
  x <- phe_pg[[2]]
  x <- (x>=thres)
  if(kick){
    s <- which(colSums(x)>0)
    s <- unique(c(1:15,s))
    x <- x[,s]
  }
  v <- colnames(x)
  m <- matrix(0,length(v),length(v),dimnames=list(v,v))
  m[1:nrow(x),1:ncol(x)] <- x
  g <- graph_from_adjacency_matrix(t(m))
  V(g)$color <- rep('blue',length(v))
  V(g)$color[1:15] <- 'red'
  plot(g,main=main,
       edge.arrow.size=.1,
       vertex.size=3,
       vertex.label.cex=.8,
       edge.width=.1)
}
plotnet2_frmomat <- function(rlt,main=""){
  g <- graph_from_adjacency_matrix(t(rlt))
  V(g)$color <- rep('blue',ncol(rlt))
  V(g)$color[1:15] <- 'red'
  plot(g,main=main,
       edge.arrow.size=.1,
       vertex.size=3,
       vertex.label.cex=.8,
       edge.width=.1)
}

################################
# Phe - Pathgroup
################################

phe <- phenet$data[,1:15]
pheadj <- phenet$adj[1:15,1:15]
phe_pg <- sparse_2sem(Y=phe,Y.fixed=pheadj,
                      X=as.matrix(do.call(cbind,pcas)),xsn=sapply(pcas,ncol),
                      lambda=0.05,times=10)
colSums(phe_pg[[2]][,-1:-15]>0)[c(18,37)]
plotnet2(phe_pg,T)


################################
# Phe - Pathgroup i
################################

length(ifrom <- which(colSums(phe_pg[[2]][,-1:-15]>=.8)>0))
ifrom <- unique(c(ifrom,18,37))
test <- lapply(ifrom,function(i){
  print(i)
  # i_sel <- unique(c(i,18,37))
  set.seed(1)
  i_sel <- i
  X.data <- do.call(c,rlt_p2pinp[[1]][i_sel])
  xsn <- sapply(X.data,ncol)
  names(xsn) <- unlist(lapply(rlt_p2pinp[[1]][i_sel],names))
  
  phe_pgi <- sparse_2sem(Y=phe,Y.fixed=pheadj,
                         X=as.matrix(do.call(cbind,X.data)),xsn=xsn,
                         lambda=ifelse(i%in%c(18,37),0.03,0.04),times=10, 
                         stability = ifelse(i%in%c(18,37),0.7,0.8))
  plotnet2(phe_pgi,T,ifelse(i%in%c(18,37),0.7,0.8),names(rlt_p2pinp[[1]])[i])
  return(phe_pgi[[3]])
})

##########################
# i_sel <- c(18,37)
# X.data <- do.call(c,rlt_p2pinp[[1]][i_sel])
# xsn <- sapply(X.data,ncol)
# names(xsn) <- unlist(lapply(rlt_p2pinp[[1]][i_sel],names))
# 
# phe_pgi <- sparse_2sem(Y=phe,Y.fixed=pheadj,
#                        X=as.matrix(do.call(cbind,X.data)),xsn=xsn,
#                        lambda=0.02,times=10)
# plotnet2(phe_pgi,T,.8,names(rlt_p2pinp[[1]])[i])
# test[[length(test)+1]] <- phe_pgi[[3]]
##################

test2 <- unique(do.call(rbind,lapply(test,function(x)x))[,c(2,4)])
v <- c(colnames(phe),unlist(lapply(rlt_p2pinp[[1]],names)))
rlt <- matrix(0,length(v),length(v),dimnames=list(v,v))
for(i in 1:nrow(test2)){
  rlt[match(test2[i,1],v),match(test2[i,2],v)] <- 1
}
plotnet2_frmomat(rlt[which(colSums(rlt)>0|rowSums(rlt>0)),which(colSums(rlt)>0|rowSums(rlt>0))])
dag_pheno_path <- rlt
save(dag_pheno_path,file='rlt_p2phe.rda')

