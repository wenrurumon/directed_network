rm(list=ls())
library(reshape)
library(dplyr)
library(igraph)

mat2df <- function(mat){
  rlt <- filter(melt(mat),value>0)[,1:2,drop=F]
  colnames(rlt) <- c('to','from')
  rlt
}
themodel <- function(iy,ix,data){
  data <- data.matrix(data)
  y <- data[,iy]
  if(length(ix)==0){
    res <- y
  } else {
    x <- data[,ix,drop=F]
    res <- lm(y~x)$residual
  }
  res <- cbind(res)
  colnames(res) <- colnames(data)[iy]
  res
}

#methylation to methylation in clusters
load("C:/Users/zhu2/Documents/getpathway/model20170215/methylation_net/methylation_site_network_new.rda")
rlt.m2m <- list(data=data.grp,adj=do.call(c,rlt2))
rlt.m2m$adj <- rlt.m2m$adj[sapply(rlt.m2m$adj,is.matrix)]
rlt.m2m$df <- do.call(rbind,lapply(rlt.m2m$adj,mat2df))
rlt.m2m$df <- cbind(to=paste0('m',rlt.m2m$df[,1]),from=paste0('m',rlt.m2m$df[,2]))
mdata <- do.call(cbind,do.call(c,data.grp))
colnames(mdata) <- paste0('m',colnames(mdata))

gene <- unique(colnames(mdata))
df <- rlt.m2m$df
mdata.res <- lapply(gene,function(genex){
  iy <- which(gene==genex)
  print(iy)
  geney <- df[which(df[,1]==genex),2]
  ix <- which(gene%in%geney)
  themodel(iy,ix,data=mdata)
})
save(mdata.res,file='C:/Users/zhu2/Documents/getpathway/model20170215/methylation_net/mdata_res.rda')
