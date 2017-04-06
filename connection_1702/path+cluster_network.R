
library(Rcpp)
rm(list=ls())
load("C:/Users/zhu2/Documents/getpathway/model20170215/rlt_p2pcp.rda")
load("C:/Users/zhu2/Documents/getpathway/model20170215/rlt_p2pinp.rda")
rlt_p2pcp <- rlt_p2pcp[[2]]
rlt_p2pinp <- rlt_p2pinp[[2]]

load("C:/Users/zhu2/Documents/getpathway/model20170215/summary/rlt_p2pfinal.rda")
rlt_p2pfinal <- rlt_p2pfinal[[2]]

setwd('C:\\Users\\zhu2\\Documents\\getpathway\\model20170215\\expression_clustering')
source('C:\\Users\\zhu2\\Documents\\sample_grouplasso_network\\scr\\grpsem.R')
load('expr_cluster.rda')
load('undnet633_4.rda')

setwd('C:\\Users\\zhu2\\Documents\\sample_grouplasso_network')
sourceCpp("scr\\score_function_regression.cpp")
sourceCpp("scr\\simple_cycle.cpp")
sourceCpp("scr\\initial_sem.cpp")
source("scr\\local_cnif_macro.R")
source("scr\\grpsem.R")

qpca2 <- function(x){
  x <- qpca(x,rank=which(qpca(x)$prop>=0.9)[1])
  x$X[,1:which(x$prop>=0.9)[1],drop=F]
}

#############################

Y <- expr_cluster
Ym <- sapply(expr_cluster,function(x){x[,1]})

Yhc1 <- lapply(lapply(rlt_p2pinp,colnames),function(x){
  which(names(expr_cluster)%in%x)
})
Yhc1 <- unlist(lapply(1:44,function(i){
  x <- rep(i,length(Yhc1[[i]]))
  names(x) <- Yhc1[[i]]
  x
}))
Yhc1 <- Yhc1[order(as.numeric(names(Yhc1)))]

Yhc2 <- hclust(dist(t(Ym[,-1:-299])))
Yhc2 <- cutree(Yhc2,15)+max(Yhc1)
Yhc <- c(Yhc1,Yhc2)
names(Yhc) <- names(Y)

Y2 <- sapply(1:length(unique(Yhc)),function(i){
  print(i)
  x <- do.call(cbind,Y[Yhc==i])
  x2 <- try(qpca(x,rank=which(qpca(x)$prop>=0.9)[1]))
  if(!is.list(x)){
    x<-qpca(x)
  }else{
    x <- qpca(x,rank=which(qpca(x)$prop>=0.9)[1])
  }
  # x$X[,1:which(x$prop>=0.9)[1],drop=F]
  x$X[,1]
})
# names(Y2) <- paste0('mcluster',1:length(Y2))
colnames(Y2) <- paste0('mcluster',1:ncol(Y2))

#############################

Y2sem <- sparse_2sem(Y=Y2,lambda=0.25)
# Y2sem <- grpsem(Y=Y2,lambda=0.5)
plot(igraph::graph_from_adjacency_matrix(Y2sem[[1]]),
     vertex.size=4,
     edge.arrow.size=.2)
Y2sem[[1]][1:44,1:44] <- rlt_p2pcp
plot(igraph::graph_from_adjacency_matrix(Y2sem[[1]]),
     vertex.size=4,
     edge.arrow.size=.2)

#############################
setwd('C:/Users/zhu2/Documents/hug')
source('codes/sparse_2sem_final.R')
source('codes/local_cnif_macro.R')
source('codes/flm_and_cca.R')
source('codes/CNIF_grouplasso.R')
source('codes/CNIF.R')
sourceCpp("codes/score_function_regression.cpp")
sourceCpp("codes/simple_cycle.cpp")
sourceCpp("codes/initial_sem.cpp")

gc()
Y2cnif <- CNIF(data=Y2,init.adj=Y2sem[[1]],max_parent=3)
plot(igraph::graph_from_adjacency_matrix(Y2cnif),
     vertex.size=4,
     edge.arrow.size=.2)
Y2cnif[1:44,1:44] <- rlt_p2pcp
is.dag(graph_from_adjacency_matrix(Y2cnif))
plot(igraph::graph_from_adjacency_matrix(t(Y2cnif)),
     vertex.size=4,
     edge.arrow.size=.2)

##########################

group_net <- lapply(1:59,function(i){
  print(i)
  Yi <- Ym[,which(Yhc==i),drop=F]
  Xi <- Ym[,Yhc%in%which(Y2cnif[i,]>0),drop=F]
  if(ncol(Yi)==1){
    out <- rbind(c(0,slim(Y=Yi,X=Xi,lambda=0.05)$beta!=0))
    dimnames(out) <- list(colnames(Yi),c(colnames(Yi),colnames(Xi)))
    YXi <- list(out)
  } else if(ncol(Xi)==0){
    YXi <- sparse_2sem(Y=Yi,lambda=0.05)
  } else {
    YXi <- sparse_2sem(Y=Yi,X=Xi,lambda=0.05)  
  }
  YXi[[1]]
})
fnet <- matrix(0,663,663,dimnames=list(names(Y),names(Y)))
for(i in 1:59){
  fneti <- group_net[[i]]
  fnet[match(rownames(fneti),rownames(fnet)),match(colnames(fneti),colnames(fnet))] <- fneti
}
fnet[match(rownames(rlt_p2pfinal),rownames(fnet)),match(colnames(rlt_p2pfinal),colnames(fnet))] <- rlt_p2pfinal

fnet2 <- fnet
dimnames(fnet2) <- list(substr(colnames(fnet),1,regexpr(' ',paste0(colnames(fnet),' '))),substr(colnames(fnet),1,regexpr(' ',paste0(colnames(fnet),' '))))

plot(igraph::graph_from_adjacency_matrix(t(fnet2)),
     vertex.size=4,
     edge.arrow.size=.2)
