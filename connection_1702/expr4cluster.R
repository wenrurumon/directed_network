
rm(list=ls())
library(reshape)
library(dplyr)
library(igraph)
library(flare)

slimsel <- function(y,x,lambda,rho=1,method='lasso',verbose=F,names=T){
  out <- slim(X=x,Y=y,lambda=lambda,rho=rho,method='lasso',verbose=F)
  sel <- which(out$beta>0)
  if(names){
    sel <- colnames(x)[sel]
  }
  return(sel)
}

#get summary
#https://github.com/wenrurumon/directed_network/blob/master/summary/pathway2phe.summary
s1 <- read.csv('clipboard')
#https://github.com/wenrurumon/directed_network/blob/master/summary/pathway2disease.summary
s2 <- readLines('clipboard')
s2 <- do.call(rbind,strsplit(s2,'\t'))[-1,]

#getdata
load('C:/Users/zhu2/Documents/getpathway/model20170215/rlt_cluster_network.rda')
rlt.g2g <- list(data=exprincluster,adj=rlt_all_split,adj_c2c=cluster_network)
load("C:/Users/zhu2/Documents/getpathway/model20170215/sample_methy_snp/disease.rda")
load("C:/Users/zhu2/Documents/getpathway/model20170215/summary/D6.rda")
disease <- D6[match(rownames(disease),rownames(D6)),]
colnames(disease) <- unique(paste(s2[,2]))

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
gene2phe <- do.call(rbind,lapply(phe,function(pi){
    y <- data.matrix(phenet$res[,which(colnames(phenet$res)==pi),drop=F])
    x <- rlt.g2g$data[which(names(rlt.g2g$data)%in%paste(filter(s1,to==pi)$from))]
    cbind(from=unlist(lapply(x,slimsel,y=y,lambda=0.05)),to=pi)
  })
)
#connect expression to disease
d <- unique(paste(s2[,2]))
gene2d <- do.call(rbind,lapply(d,function(pi){
  y <- data.matrix(disease[,which(colnames(disease)==pi),drop=F])
  x <- rlt.g2g$data[which(names(rlt.g2g$data)%in%s2[s2[,2]==pi,1])]
  cbind(from=unlist(lapply(x,slimsel,y=y,lambda=0.05)),to=pi)
}))


