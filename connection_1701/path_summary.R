
rm(list=ls())
library(igraph)
setwd('C:\\Users\\zhu2\\Documents\\getpathway\\model20170126')
for(i in dir(pattern='rda')){load(i)}

plotnet <- function(x,mode='directed',main=''){
  diag(x) <- 0
  plot(graph_from_adjacency_matrix(t(x),mode=mode),
       edge.arrow.size=.3,
       vertex.size=3,
       vertex.label.cex=1,
       edge.width=.1,
       main=main)
}

# par(mfrow=c(1,2))
# for(i in grep('Signal',names(rlt_p2pfinal[[1]]))){
#   temp <- rlt_p2pfinal[[1]][[i]]
#   plotnet(temp[,1:nrow(temp)])
# }
# 
# g <- graph_from_adjacency_matrix(t(temp[,1:nrow(temp)]))
# lapply(1:5,function(i){all_simple_paths(g,from=V(g)[i])})
# plotnet(temp[,1:nrow(temp)])
# 
# tmp <- matrix(0,0,2)
# for(i in 1:nrow(rlt_p2pfinal[[2]])){
#   for(j in 1:ncol(rlt_p2pfinal[[2]])){
#     if(rlt_p2pfinal[[2]][i,j]>0){
#       tmp <- rbind(tmp,c(rownames(rlt_p2pfinal[[2]])[i],colnames(rlt_p2pfinal[[2]])[j]))
#     }    
#   }
# }
# colnames(tmp) <- c('y','x')
# write.csv(tmp,'p2pfinal.csv',row.names=F)

plotnet(rlt_p2pfinal[[2]])
par(mfrow=c(1,1))
for(i in grep('Signal',names(rlt_p2pfinal[[1]]))){
  temp <- rlt_p2pfinal[[1]][[i]]
  plotnet(temp[,1:nrow(temp),drop=F],main=names(rlt_p2pfinal[[1]])[i])
}
par(mfrow=c(2,2))
for(i in 1:length(rlt_p2pfinal[[1]])){
  print(i)
  temp <- rlt_p2pfinal[[1]][[i]]
  plotnet(temp[,1:nrow(temp),drop=F],main=names(rlt_p2pfinal[[1]])[[i]])
}

##############################
load('pathwaymap.rda')
temp <- t(rlt_p2pfinal[[2]])
g <- graph_from_adjacency_matrix(temp)

temp <- rownames(rlt_p2pfinal[[1]][[18]])
pathwaymap <- pathwaymap[names(pathwaymap)%in%temp]
pathwaymap <- lapply(pathwaymap, function(x){unique(x[x%in%temp])})

temp <- lapply(1:length(pathwaymap),function(i){
  print(i)
  pathi <- unique(pathwaymap[[i]])[!grepl(names(pathwaymap)[i],unique(pathwaymap[[i]]))]
  node.from <- V(g)[grep(names(pathwaymap)[i],names(V(g)))]
  node.to <- V(g)[grep(paste(pathwaymap[[i]],collapse="|"),names(V(g)))]
  paths <- all_shortest_paths(g,node.from,node.to)$res
  list(paths=paths,ends=unique(sapply(paths,function(x){x[length(x)]})))
})
sum(sapply(temp,function(x){length(x$ends)}))/sum(sapply(pathwaymap,length))


