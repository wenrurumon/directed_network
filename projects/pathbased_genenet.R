
rm(list=ls())
setwd('C:\\Users\\zhu2\\Documents\\network_final')
# load("~/getpathway/gene39761.rda.RData")
# Y.exp <- expinpath[-3]
rm(list=ls()[-grep('Y.exp',ls())])
load('pathnet.rda'); load('rlt_expnetwork.rda')

#inputdata
setwd('C:\\Users\\zhu2\\Documents\\signaling\\codes\\')
source('sparse_2sem_final.R')
source('local_cnif_macro.R')
source('CNIF.R')
sourceCpp("score_function_regression.cpp")
sourceCpp("simple_cycle.cpp")
sourceCpp("initial_sem.cpp")
setwd('C:\\Users\\zhu2\\Documents\\getpathway')
library(flare)
library(grplasso)
library(data.table)
library(dplyr)
library(igraph)

#
i <- 1
pathi <- names(rlt)[i]
x.path <- colnames(pathnet$pathnet)[pathnet$pathnet[grep(pathi,rownames(pathnet$pathnet)),,drop=F]>=.1]
Y <- as.matrix(rlt[[i]]$data)
Y.fixed <- rlt[[i]]$adj$eq_matrix
Xs <- lapply(match(x.path,names(rlt)),function(i){
  X <- rlt[[i]]$data
  as.matrix(X)
})

temps <- lapply(Xs,function(x){
  print('newmodel')
  temp <- sparse_2sem(Y=Y,Y.fixed=Y.fixed,X=x,lambda=0.2)
  temp$eq_matrix
})
