rm(list=ls())
library(reshape)
library(dplyr)
library(igraph)

mat2df <- function(mat){
  rlt <- filter(melt(mat),value>0)[,1:2,drop=F]
  colnames(rlt) <- c('to','from')
  rlt
}
sumnet <- function(x){
  x <- unique(x)
  v <- unique(as.vector(x))
  c(v = length(v), e = nrow(x))
}
togene <- function(x){
  type <- substr(x,1,1)
  x <- strsplit(x,':')[[1]]
  
  
}

#methylation to methylation in clusters
load("C:/Users/zhu2/Documents/getpathway/model20170215/methylation_net/methylation_site_network_new.rda")
rlt.m2m <- list(data=data.grp,adj=do.call(c,rlt2))
rlt.m2m$adj <- rlt.m2m$adj[sapply(rlt.m2m$adj,is.matrix)]
rlt.m2m$df <- do.call(rbind,lapply(rlt.m2m$adj,mat2df))
rlt.m2m$df <- cbind(to=paste0('m',rlt.m2m$df[,1]),from=paste0('m',rlt.m2m$df[,2]))

#expression to expression cross all
load('C:/Users/zhu2/Documents/getpathway/model20170215/rlt_cluster_network.rda')
rlt.g2g <- list(data=exprincluster,adj=rlt_all_split,adj_c2c=cluster_network)
rlt.g2g$df <- (do.call(rbind,lapply(rlt.g2g$adj,mat2df)))
rlt.g2g$df <- cbind(to = paste(rlt.g2g$df[,1]), from = paste( rlt.g2g$df[,2]))

#qtl, methylation or snp to expression
load('C:\\Users\\zhu2\\Documents\\getpathway\\model20170215\\tacc\\qtlresult.rda')
rlt.qtl <- list(rlt=rlt)
rlt.qtl$df <- rbind(
  cbind(to=paste(rlt[[1]]$X1),from=paste0('s',rlt[[1]]$X2)),
  cbind(paste(rlt[[2]]$X1),paste0('m',rlt[[2]]$X2)),
  cbind(paste(rlt[[3]]$X2),paste(rlt[[3]]$X1))
)

sumnet(rlt.m2m$df)
sumnet(rlt.g2g$df)
sumnet(rlt.qtl$df)
459248/(53011*53011)

df.m2m <- rlt.m2m$df
df.m2m[,1] <- paste(substr(df.m2m[,1],1,1),sapply(strsplit(df.m2m[,1],':'),function(x){x[[2]]}),sep=':')
df.m2m[,2] <- paste(substr(df.m2m[,2],1,1),sapply(strsplit(df.m2m[,2],':'),function(x){x[[2]]}),sep=':')
df.g2g <- rlt.g2g$df

