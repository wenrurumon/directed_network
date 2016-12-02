
library(igraph)
library(dplyr)
rm(list=ls())
setwd('C:\\Users\\zhu2\\Documents\\getpathway')
raw <- as.matrix(read.csv('20161026network.csv'))
# g <- raw[grepl('pheno::',raw[,1])&grepl('pheno::',raw[,2]),]
g <- cbind(raw[,2],raw[,1])
g <- graph_from_data_frame(g)
plot(g)
#find path
from = "geno::CREBBP"
to = 'disease::AD'
getpath <- function(g,from,to){
      key <- colnames(as_adjacency_matrix(g))
      from <- grep(from,key)
      to <- grep(to,key)
      igraph::all_simple_paths(g,from,to)
}
rlt <- getpath(g,from,to)
