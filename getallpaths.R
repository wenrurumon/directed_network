
library(igraph)
library(dplyr)

rm(list=ls())

setwd('C:\\Users\\zhu2\\Documents\\getpathway')
raw <- as.matrix(read.csv('20161026network.csv'))
g <- cbind(raw[,2],raw[,1])
g <- graph_from_data_frame(g)
plot(g)

#find path
getpath <- function(g,from,to){
      key <- colnames(as_adjacency_matrix(g))
      from <- grep(from,key)
      to <- grep(to,key)
      out <- igraph::all_simple_paths(g,from,to)
      out <- lapply(out,function(x) paste(names(x),collapse=' -> '))
      unlist(out)
}

rlt <- lapply(c('geno::ELF3','geno::CREBBP','geno::POU3F2'),function(x){
  print(x)
  getpath(g,x,'disease::AD')
})

key <- colnames(as_adjacency_matrix(g))
geno <- grep('geno::',key,value=T)
d <- grep('disease::',key,value=T)

test <- sapply(geno,function(from){
  sapply(d,function(to){
    length(getpath(g,from,to))
  })
})

test2 <- getpath(g,'geno::CREBBP','disease::AD')
