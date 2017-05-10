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

#methylation to methylation in clusters
load("C:/Users/zhu2/Documents/getpathway/model20170215/methylation_net/methylation_site_network_new.rda")
rlt.m2m <- list(data=data.grp,adj=do.call(c,rlt2))
rlt.m2m$adj <- rlt.m2m$adj[sapply(rlt.m2m$adj,is.matrix)]
rlt.m2m$df <- do.call(rbind,lapply(rlt.m2m$adj,mat2df))
rlt.m2m$df <- cbind(to=paste0('m',rlt.m2m$df[,1]),from=paste0('m',rlt.m2m$df[,2]))
rdf.m2m <- rlt.m2m$df

#expression to expression cross all
load('C:/Users/zhu2/Documents/getpathway/model20170215/rlt_cluster_network.rda')
rlt.g2g <- list(data=exprincluster,adj=rlt_all_split,adj_c2c=cluster_network)
rlt.g2g$df <- (do.call(rbind,lapply(rlt.g2g$adj,mat2df)))
rlt.g2g$df <- cbind(to = paste(rlt.g2g$df[,1]), from = paste( rlt.g2g$df[,2]))
rdf.g2g <- rlt.g2g$df

#qtl, methylation or snp to expression
setwd('C:\\Users\\zhu2\\Documents\\getpathway\\model20170215\\tacc\\')
for(i in dir(pattern='df')){load(i)}
rdf.m2g <- gsub('::',':',df_m2eres)
rdf.s2g <- gsub('::',':',df_s2eres)
rdf.s2m <- gsub('::',':',df_s2mres)

#clipboard from github
rdf.g2d <- read.table('clipboard',header=T,as.is=c('from','to'))
# rdf.g2d <- cbind(rdf.g2d[,2],rdf.g2d[,1])
rdf.g2p <- read.csv('clipboard',header=T,as.is=c('from','to'))
rdf.g2p <- cbind(rdf.g2p[,2],rdf.g2p[,1])
rdf.p2d <- read.table('clipboard',header=T,as.is=c('from','to'))
rdf.p2d <- cbind(rdf.p2d[,2],rdf.p2d[,1])
rdf.p2p <- read.table('clipboard',header=T,as.is=c('from','to'))
rdf.p2p <- cbind(rdf.p2p[,2],rdf.p2p[,1])
rdf.lr <- read.table('clipboard',header=T,as.is=c('from','to'))
rdf.lr <- cbind(rdf.lr[,2],rdf.lr[,1])

#RESULT
rlt <- ls(pattern='rdf')
rlt <- lapply(rlt,function(x){
  x <- eval(parse(text=x))
  colnames(x) <- c('to','from')
  x
})
names(rlt) <- ls(pattern='rdf')
setwd('C:\\Users\\zhu2\\Documents\\getpathway\\model20170215\\summary_201705\\')
load('rdf.rda')
# save(rlt,file='rdf.rda')

#Summary
rlt2 <- as.matrix(do.call(rbind,rlt))
rlt2 <- rlt2[rlt2[,1]!=rlt2[,2],]
x <- rlt2
sel <- substr(x,1,1)%in%c('m','s')&grepl(":",substr(x,3,4))
xsel <- x[sel]
xsel <- paste0(substr(xsel,1,1),substr(xsel,regexpr(':',xsel),nchar(xsel)))
rlt2 <- as.matrix(rlt2)
rlt2[sel] <- xsel
rlt2 <- unique(rlt2)

g <- igraph::graph_from_data_frame(data.frame(rlt2[,2],rlt2[,1]))
shortest_paths(g,from='POU3F2',to='AD')
all_simple_paths(g,from='m:POU3F2')


