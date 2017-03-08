g <- graph_from_adjacency_matrix(t(xmat),'undirected')
test <- fastgreedy.community(g)
test <- table(test$membership,xnode[,2])
test2 <- matrix(0,nrow(test),ncol(test),dimnames=dimnames(test))
for(i in 1:nrow(test2)){
  for(j in 1:ncol(test2)){
    m <- colSums(test)[j]
    n <- sum(test) - m
    k <- rowSums(test)[i]
    x <- test[i,j]
    test2[i,j] <- dhyper(x,m,n,k,log=F)
  }
}
