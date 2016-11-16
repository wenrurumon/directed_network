qpca <- function(A,rank=0,ifscale=TRUE){
  if(ifscale){A <- scale(as.matrix(A))[,]}
  A.svd <- svd(A)
  if(rank==0){
    d <- A.svd$d
  } else {
    d <- A.svd$d-A.svd$d[min(rank+1,nrow(A),ncol(A))]
  }
  d <- d[d > 1e-8]
  r <- length(d)
  prop <- d^2; info <- sum(prop)/sum(A.svd$d^2);prop <- cumsum(prop/sum(prop))
  d <- diag(d,length(d),length(d))
  u <- A.svd$u[,1:r,drop=F]
  v <- A.svd$v[,1:r,drop=F]
  x <- u%*%sqrt(d)
  y <- sqrt(d)%*%t(v)
  z <- x %*% y
  rlt <- list(rank=r,X=x,Y=y,Z=x%*%y,prop=prop,info=info)
  return(rlt)
}
pca <- function(X){
  X <- scale(as.matrix(X))
  m = nrow(X)
  n = ncol(X)
  Xeigen <- svd(X)
  value <- (Xeigen$d)^2/m
  value <- cumsum(value/sum(value))
  score <- X %*% Xeigen$v
  mat <- Xeigen$v
  list(score=score,prop=value,mat=mat)
}

plotnet <- function(x,mode='undirected'){
  diag(x) <- 0
  plot(graph_from_adjacency_matrix(t(x),mode=mode),
       edge.arrow.size=.1,
       vertex.size=3,
       vertex.label.cex=1,
       edge.width=.1)
}
fc <- function(x){
  w<-as.vector(t(x))[t(x)>0]
  x <- graph_from_adjacency_matrix(x>0,mode='undirected')
  fc <- membership(fastgreedy.community(x,weight=w))
  fc[] <- match(fc,unique(fc))
  fc
}
plotclust <- function(x,membership=NULL){
  G <- graph_from_adjacency_matrix(x>0)
  if(is.null(membership)){membership=rep(1,ncol(x))}
  plot(create.communities(G, membership), 
       # as.undirected(G), 
       as.directed(G),
       layout=layout.kamada.kawai(as.undirected(G)),
       edge.arrow.size=.2,
       vertex.size=3,
       vertex.label.cex=1,
       edge.width=.2)
}
