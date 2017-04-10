
rm(list=ls())

library(Rcpp)
library(grplasso)
library(data.table)
library(dplyr)
library(igraph)
library(RcppArmadillo)
library(hash)
library(rstackdeque)
library(Rsymphony)

setwd('C:\\Users\\zhu2\\Documents\\sample_grouplasso_network')
sourceCpp("scr\\score_function_regression.cpp")
sourceCpp("scr\\simple_cycle.cpp")
sourceCpp("scr\\initial_sem.cpp")
source("scr\\local_cnif_macro.R")
source("scr\\grpsem.R")

load('expressiondata.rda')

par(mfrow=c(1,2))

Y1 <- Y.expression[1:5]
system.time(test1 <- model(Y1,x.graph=NULL,prop=0.8,lambda=0.5,max.parent = 3))
plot(graph_from_adjacency_matrix(t(test1)))

Y2 <- lapply(Y1,function(x) scale(x[,1:30,drop=F]))
system.time(test2 <- model(Y2,x.graph=NULL,prop=0.8,lambda=0.5,max.parent = 3))
plot(graph_from_adjacency_matrix(t(test2)))
