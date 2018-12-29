rm(list=ls())
library(igraph)
library(WGCNA)
library(grplasso)
setwd("/Users/wenrurumon/Downloads/huzixin/")
dir()
load('expr_for_WGCNA.rda')

plotnet <- function(x,mode='directed',main=''){
  diag(x) <- 0
  plot(graph_from_adjacency_matrix(t(x),mode=mode),
       edge.arrow.size=.3,
       vertex.size=1,
       vertex.label.cex=1,
       edge.width=.1,
       main=main)
}

x <- expr_all
TOM <- TOMsimilarityFromExpr(x,power=6)
dimnames(TOM) <- list(colnames(x),colnames(x))
diag(TOM) <- 0
# cyt <- exportNetworkToCytoscape(TOM,nodeNames=names(datExpr),threshold=0.02)
cenet <- (TOM>=quantile(TOM,0.95));sum(cenet)
par(mfrow=c(1,1)); plotnet(cenet>0,'undirected')

set.seed(1234);test <- sample(1:447,47)
expr_train <- expr_all[-test,]
expr_test <- expr_all[test,]

pred_tom_test <- sapply(1:ncol(expr_train),function(i){
  print(i)
  yi_train <- expr_train[,i]
  xi_train <- expr_train[,-i]
  if(ncol(xi_train)==0){
    return(rep(mean(yi_train),nrow(data_test)))
  } else {
    xi_test <- expr_test[,-i]
    xi_lm <- coef(lm(yi_train~xi_train))
    as.numeric(data.matrix(cbind(1,xi_test)) %*% matrix(xi_lm,ncol=1))
  }
})

tom_rlt <- cor(pred_tom_test,expr_test)
diag(tom_rlt)

tom_res <- (abs(pred_tom_test-expr_test)) 
mean(colMeans(tom_res))
quantile(colMeans(tom_res))
