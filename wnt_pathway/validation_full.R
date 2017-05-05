# rm(list=ls())
setwd('C:\\Users\\zhu2\\Documents\\signaling\\codes\\')
source('sparse_2sem_final.R')
source('local_cnif_macro.R')
source('CNIF.R')
sourceCpp("score_function_regression.cpp")
sourceCpp("simple_cycle.cpp")
sourceCpp("initial_sem.cpp")
setwd('C:\\Users\\zhu2\\Documents\\getpathway')
load("~/getpathway/gene39761.rda.RData")
load('C:/Users/zhu2/Documents/getpathway/model20170215/rlt_p2pinp.rda')
library(flare)
library(grplasso)
library(data.table)
library(dplyr)
library(igraph)
# rlt.jr <- as.matrix(read.table('clipboard',header=T))

############################
# SEM
############################

expr <- expinpath$`Wnt signaling pathway`
# expr <- expr[,colnames(expr)%in%as.vector(rlt.jr)]
gc();sem1 <- sparse_2sem(expr,times=100,lambda=1)
sum(sem1[[1]]>=0.8)
gc();dag <- CNIF(data=expr,init.adj=(sem1[[1]]>=.8),max_parent=3)
sem <- sparse_2sem(expr,Y.fixed=dag)[[2]][,c(2,4)]

############################
# Validation
############################

# ref <- read.table('clipboard',header=T)
g.ref <- graph_from_data_frame(ref)
from <- sem[1,1]
to <- sem[1,2]

valii <- function(from,to,g.ref){
  ifrom <- which(names(V(g.ref))==from)
  ito <- which(names(V(g.ref))==to)
  shortest_paths(g.ref,V(g.ref)[ifrom],V(g.ref)[ito])
}
test <- lapply(1:nrow(sem),function(i){
  from <- sem[i,1]
  to <- sem[i,2]
  valii(from,to,g.ref)$vpath[[1]]
})
mean(sapply(test,length)>0)
