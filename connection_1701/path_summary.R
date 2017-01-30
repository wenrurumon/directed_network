
rm(list=ls())
setwd('C:\\Users\\zhu2\\Documents\\getpathway\\model20170126')
for(i in dir(pattern='rda')){load(i)}

plotnet <- function(x,mode='directed'){
  diag(x) <- 0
  plot(graph_from_adjacency_matrix(t(x),mode=mode),
       edge.arrow.size=.3,
       vertex.size=3,
       vertex.label.cex=1,
       edge.width=.1)
}

plotnet(rlt_p2pfinal[[2]])
