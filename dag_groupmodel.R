
rm(list=ls())
gc()

# load('network_per_group.rda')
setwd('C:\\Users\\zhu2\\Documents\\signaling\\codes\\')
source('sparse_2sem_final.R')
source('local_cnif_macro.R')
source('CNIF.R')
sourceCpp("score_function_regression.cpp")
sourceCpp("simple_cycle.cpp")
sourceCpp("initial_sem.cpp")
setwd('C:\\Users\\zhu2\\Documents\\getpathway')
load("~/getpathway/gene39761.rda.RData")
load("C:/Users/zhu2/Documents/getpathway/ptwmap2.rda")
library(grplasso)
library(data.table)
library(dplyr)
library(igraph)

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

qexpinpath <- lapply(expinpath,function(x){
  if(ncol(x)==1){
    return(x)
  }else{
    x <- qpca(x,which(pca(x)$prop>=0.9)[1])
    x$X[,1:which(x$prop>=0.9)[1],drop=F]
  }  
}
)
plotnet <- function(x,mode='directed'){
  diag(x) <- 0
  plot(graph_from_adjacency_matrix(t(x),mode=mode),
       edge.arrow.size=.1,
       vertex.size=3,
       vertex.label.cex=1,
       edge.width=.1)
}
pcainpath <- lapply(qexpinpath,function(x){pca(x)})
for(i in 1:length(pcainpath)){
  colnames(pcainpath[[i]]$score) <- paste(ptwmap[match(names(pcainpath),ptwmap[,2]),1][i],1:ncol(pcainpath[[i]]$score),sep="_")
}
pcas <- lapply(pcainpath,function(x){x[[1]]})[-3]

#################################################
# sub Macros
#################################################

equationj <- function(j,x.input,lambda=0.5){
  y <- x.input[[j]]
  x <- x.input[-j]
  #process data
  processx <- function(x,y){
    temp <- matrix(0,nrow=nrow(y)*ncol(y),ncol=ncol(x)*ncol(y))
    for(i2 in 1:ncol(y)){
      temp[(i2-1)*nrow(y)+1:nrow(y),(i2-1)*ncol(x)+1:ncol(x)] <- scale(x)
    }
    temp
  }
  x.p <- lapply(x.input[-j],function(x){
    processx(x,x.input[[j]])
  })
  for(i2 in 1:length(x.p)){colnames(x.p[[i2]]) <- paste(names(x.p)[i2],1:ncol(x.p[[i2]]),sep='_')}
  y <- cbind(scale(as.vector(y)))
  #group lasso model
  lambda <- lambdamax(x=cbind(1,do.call(cbind,x.p)),y=y, 
                      index=c(NA,rep(1:length(x),sapply(x.p,ncol))), 
                      penscale = sqrt, model = LinReg(),
                      center=TRUE,standardized=TRUE) * 0.5^(1/lambda-1)
  fit <- grplasso(x=cbind(1,do.call(cbind,x.p)),y=y,
                  index=c(NA,rep(1:length(x),sapply(x.p,ncol))),lambda=lambda,model=LinReg(),
                  penscale = sqrt,
                  control = grpl.control(update.hess = "lambda", trace = 0))
  temp <- rep(0,length(x.input))
  temp[-j] <- tapply(coef(fit)[-1]!=0,rep(1:length(x),sapply(x.p,ncol)),function(x){any(x!=0)})
  temp
}
subpath <- function(pathi, max.parent = 3){
  n.parent <- min(sum(pathi),max.parent)
  path0 <- rep(0,length(pathi))
  if(sum(pathi)==0){
    matrix(path0,nrow=1)
  } else {
    paths <- do.call(rbind,
                     lapply(1:n.parent,function(np){
                       sel <- matrix(which(pathi>0)[combn(sum(pathi),np)],nrow=np)
                       t(sapply(1:ncol(sel),function(seli){
                         temp <- rep(0,length(pathi))
                         temp[sel[,seli]] <- 1
                         temp
                       }))
                     })
    )
    rbind(path0,paths)
  }
}
mape <- function(adj,x.input){
  apply(adj,1,function(p){
    yi <- scale(x.input[[p[length(p)]]])
    if(sum(p!=0)==1){
      return(mean(yi^2))
    } else {
      xi <- scale(do.call(cbind,x.input[which(p[-length(p)]>0)]))
      return(mean((lm(yi~xi)$residual)^2))
    }
  })
}
ip <- function(x.graph,x.path,x.cost){
  
  x.graph <- do.call(rbind,adj)
  x.path <- adj2
  x.cost <- cost
  
  parms <- ncol(x.graph)
  x.p <- data.frame(response=x.path[,ncol(x.path)],
                    nparent=rowSums(x.path[,2:ncol(x.path)-1,drop=F]),
                    cost=x.cost)
  x.p <- as.matrix(x.p %>% group_by(response,nparent) %>% summarise(min=min(cost)))
  
  keep.path <- rep(T,nrow(x.path))
  for(i in 1:nrow(x.p)){
    if(x.p[i,2]>1){
      x.sel <- (rowSums(x.path[,2:ncol(x.path)-1,drop=F])==x.p[i,2])&(x.path[,ncol(x.path)]==x.p[i,1])
      keep.path[x.sel] <- x.cost[x.sel]<x.p[i-1,3]
    }
  }
  
  x.path <- x.path[keep.path,]
  obj <- x.cost[keep.path]
  types <- rep('I',parms)
  
  rhs   = hash()
  sense   = hash()
  A     = hash()
  
  working.node = x.path[,parms+1]
  working.parent = x.path[,1:parms]
  for( i in 1:parms ){
    cons.name = paste("one parent set",i)
    rhs[[cons.name]] = 1
    sense[[cons.name]] = "=="
    A[[cons.name]]   = (working.node == i)
  }
  
  cons.names = keys(A)
  A = t(values(A,cons.names))
  sense = values(sense,cons.names)
  rhs = values(rhs,cons.names)
  csc.A = make.constraints(x.graph,cbind(x.path,obj))
  
  A = rbind(A, csc.A)
  sense = c(sense, rep("<=",nrow(csc.A)))
  rhs = c(rhs, rowSums(csc.A)-1)
  
  result = Rsymphony_solve_LP(obj=obj,mat=A,dir=sense,rhs=rhs,types=types,verbosity = 0 )
  dag <- x.path[result$solution==1,1:parms,drop=F]
  return(dag)
}

#################################################
# Model
#################################################

i <- 2

grpi <- unique(pathlist[,1])[i]
paths <- pathlist[pathlist[,1]==grpi,2]
input <- pcainpath[names(pcas)%in%paths]
print(list(i=i,groupname=grpi,paths=names(input)))
x.input <- lapply(input,function(x){
  x$score[,1:which(x$prop>=0.8)[1],drop=F]
})

adj <- lapply(1:length(x.input),function(j){equationj(j,x.input,lambda=0.6)})
adj2 <- do.call(rbind,lapply(1:length(adj),function(i){cbind(subpath(adj[[i]],max.parent = 3),i)}))
cost <- mape(adj2,x.input)
