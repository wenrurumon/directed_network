
library(Rcpp)
rm(list=ls())
load("C:/Users/zhu2/Documents/getpathway/model20170215/rlt_p2pinp.rda")

net <- rlt_p2pinp[[2]][[18]]
cname <- colnames(net)
cname <- substr(cname,1,regexpr(' ',cname)-1)
dimnames(net) <- list(cname,cname)

plot(igraph::graph_from_adjacency_matrix(t(net)),
     vertex.size=4,
     edge.arrow.size=.2)

net2df <- function(net){
  df <- lapply(1:nrow(net),function(i){
    responsor <- rownames(net)[i]
    predictor <- names(which(net[i,]>0))
    if(length(predictor)>0){
      return(cbind(y=responsor,x=predictor))
    } else {
      return(NULL)
    }
  })
  do.call(rbind,df)
}

write.csv(net2df(net),'C:/Users/zhu2/Documents/getpathway/model20170215/signaling_pathway.csv',row.names=F)
