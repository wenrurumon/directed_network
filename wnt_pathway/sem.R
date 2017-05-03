rm(list=ls())
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
rlt.jr <- as.matrix(read.table('clipboard',header=T))

expr <- expinpath$`Wnt signaling pathway`
expr <- expr[,colnames(expr)%in%as.vector(rlt.jr)]
gc();sem <- sparse_2sem(expr,times=100,lambda=0.6)
sum(sem[[1]]>0.8)
gc();dag <- CNIF(data=expr,init.adj=(sem[[1]]>0.8),max_parent=4)
sem <- sparse_2sem(expr,Y.fixed=dag)[[2]][,c(2,4)]

