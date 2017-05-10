rm(list=ls())
library(reshape)
library(dplyr)
library(igraph)

mat2df <- function(mat){
  rlt <- filter(melt(mat),value>0)[,1:2,drop=F]
  colnames(rlt) <- c('to','from')
  rlt
}
#expression to expression cross all
load('C:/Users/zhu2/Documents/getpathway/model20170215/rlt_cluster_network.rda')
rlt.g2g <- list(data=exprincluster,adj=rlt_all_split,adj_c2c=cluster_network)
rlt.g2g$df <- (do.call(rbind,lapply(rlt.g2g$adj,mat2df)))
rlt.g2g$df <- cbind(to = paste(rlt.g2g$df[,1]), from = paste( rlt.g2g$df[,2]))
rlt.g2g$df <- rlt.g2g$df[rlt.g2g$df[,1]!=rlt.g2g$df[,2],,drop=F]

#gene expression residual
gdata <- do.call(cbind,rlt.g2g$data)
gene <- unique(colnames(gdata))

i <- 0
expr_residual <- lapply(gene,function(geney){
  print(i<<-i+1)
  genex <- rlt.g2g$df[which(rlt.g2g$df[,1]==geney),2]
  y <- gdata[,match(geney,gene)]
  if(length(genex)==0){
    return(y)
  } else {
    x <- gdata[,match(genex,gene)]
    return(lm(y~x)$residual)
  }
})
