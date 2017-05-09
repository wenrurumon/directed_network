
rm(list=ls())
library(reshape)
library(dplyr)

mat2df <- function(mat){
  rlt <- filter(melt(mat),value>0)[,1:2,drop=F]
  colnames(rlt) <- c('to','from')
  rlt
}

#methylation to methylation in clusters
load("C:/Users/zhu2/Documents/getpathway/model20170215/methylation_net/methylation_site_network_new.rda")
rlt.m2m <- list(data=data.grp,adj=do.call(c,rlt2))
rlt.m2m$adj <- rlt.m2m$adj[sapply(rlt.m2m$adj,is.matrix)]
rlt.m2m$df <- do.call(rbind,lapply(rlt.m2m$adj,mat2df))

#expression to expression cross all
load('C:/Users/zhu2/Documents/getpathway/model20170215/rlt_cluster_network.rda')
rlt.g2g <- list(data=exprincluster,adj=rlt_all_split,adj_c2c=cluster_network)
rlt.g2g$df <- do.call(rbind,lapply(rlt.g2g$adj,mat2df))

#qtl, methylation or snp to expression
load('C:\\Users\\zhu2\\Documents\\getpathway\\model20170215\\tacc\\qtlresult.rda')
rlt.qtl <- list(rlt=rlt)
rlt.qtl$df <- rbind(
  cbind(to=paste(rlt[[1]]$X1),from=paste0('s',rlt[[1]]$X2)),
  cbind(paste(rlt[[2]]$X1),paste0('m',rlt[[2]]$X2)),
  cbind(paste(rlt[[3]]$X2),paste(rlt[[3]]$X1))
)

#merged
mrlt <- rbind(rlt.m2m$df,rlt.g2g$df,rlt.qtl$df); mrlt <- cbind(from=paste(mrlt[,2]),to=paste(mrlt[,1]))
g <- igraph::graph_from_data_frame(mrlt)
length(V(g))
