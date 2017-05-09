
rm(list=ls())
library(reshape)
library(dplyr)
library(igraph)
library(flare)

slim2 <- function(x){
  
}

#get summary
#https://github.com/wenrurumon/directed_network/blob/master/summary/pathway2phe.summary
s1 <- read.csv('clipboard')

#getdata
load('C:/Users/zhu2/Documents/getpathway/model20170215/rlt_cluster_network.rda')
rlt.g2g <- list(data=exprincluster,adj=rlt_all_split,adj_c2c=cluster_network)

#Pheno Net Residual
load('c:\\Users\\zhu2\\Documents\\getpathway\\model20170215\\phenet.rda')
colnames(phenet$adj) <- colnames(phenet$data)
phenet$res <- sapply(1:nrow(phenet$adj),function(i){
  y <- data.matrix(phenet$data[,i,drop=F])
  x <- which(phenet$adj[,i]>0)
  x <- data.matrix(phenet$data[,x,drop=F])
  if(length(x)==0){
    lm(y~1)$residual
  } else {
    lm(y~x)$residual
  }
})
colnames(phenet$res) <- colnames(phenet$data)

#connect expression to phenotype
phe <- paste(unique(s1$to))
lapply(phe,function(pi){
  
})

pi <- phe[1]

y <- data.matrix(phenet$res[,which(colnames(phenet$res)==pi),drop=F])
x <- rlt.g2g$data[which(names(rlt.g2g$data)%in%paste(filter(s1,to==pi)$from))]
lapply(x,function(xi){
  
})
