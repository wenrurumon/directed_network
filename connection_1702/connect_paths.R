
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
library(grplasso)
library(data.table)
library(dplyr)
library(igraph)

#################################################
# macros 1
#################################################

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
plotnet <- function(x,mode='directed'){
  diag(x) <- 0
  plot(graph_from_adjacency_matrix(t(x),mode=mode),
       edge.arrow.size=.3,
       vertex.size=3,
       vertex.label.cex=1,
       edge.width=.1)
}

#################################################
# macros 2
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
  gc()
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

ip_zy <- function(x.graph,x.path,x.cost){
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
ip_hzx <- function(x.graph,x.path,x.cost){
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
  csc.A = make.constraints_byigraph(x.graph,cbind(x.path,obj))
  
  if(!is.null(csc.A)){
    A = rbind(A, csc.A)
    sense = c(sense, rep("<=",nrow(csc.A)))
    rhs = c(rhs, rowSums(csc.A)-1)
    result = Rsymphony_solve_LP(obj=obj,mat=A,dir=sense,rhs=rhs,types=types,verbosity = 0 )
    dag <- x.path[result$solution==1,1:parms,drop=F]
    return(dag)
  } else {
    return(x.graph)
  }
  
}

model <- function(input,lambda=0.6 ,max.parent=3){
  gc()
  if(length(input)==1){
    dag <- matrix(0,1,1)
    dimnames(dag) <- list(names(input),names(input))
    return(list(dag=dag))
  }
  
  x.input <- lapply(input,function(x){
    # x$score[,1:which(x$prop>=0.8)[1],drop=F]
    x$score
  })
  
  adj <- lapply(1:length(x.input),function(j){equationj(j,x.input,lambda=lambda)})
  adj2 <- do.call(rbind,lapply(1:length(adj),function(i){cbind(subpath(adj[[i]],max.parent = max.parent),i)}))
  cost <- mape(adj2,x.input)
  dag <-   ip(do.call(rbind,adj),adj2,cost)
  dimnames(dag) <- list(names(x.input),names(x.input))
  # plotnet(dag)
  list(dag=dag,coef=c(lambda,max.parent))
}


model2 <- function(input,lambda,max.parent){
  
  gc()
  if(length(input)==1){
    dag <- matrix(0,1,1)
    dimnames(dag) <- list(names(input),names(input))
    return(list(dag=dag))
  }
  
  x.input <- lapply(input,function(x){
    # x$score[,1:which(x$prop>=0.8)[1],drop=F]
    x$score
  })
  
  adj <- lapply(1:length(x.input),function(j){equationj(j,x.input,lambda=lambda)})
  adj2 <- do.call(rbind,lapply(1:length(adj),function(i){cbind(subpath(adj[[i]],max.parent = max.parent),i)}))
  cost <- mape(adj2,x.input)
  dag <-   ip_zy(do.call(rbind,adj),adj2,cost)
  
  adj <- lapply(1:nrow(dag),function(i){dag[i,]})
  adj2 <- do.call(rbind,lapply(1:length(adj),function(i){cbind(subpath(adj[[i]],max.parent = max.parent),i)}))
  cost <- mape(adj2,x.input)
  dag <-   ip_hzx(do.call(rbind,adj),adj2,cost)
  
  dimnames(dag) <- list(names(x.input),names(x.input))
  plotnet(dag)
  
  return(dag)  
}

get_all_circles <- function(mat,score){
  g <- graph_from_adjacency_matrix(mat)
  do.call(c,lapply(1:ncol(mat),function(nodei){
    # print(nodei/ncol(mat))
    lapply(all_simple_paths(g,nodei,unique(score[score[,nodei]==1,ncol(mat)+1])),function(x){c(x,x[1])})
  }))
}
pastes <- function(x){
  xi <- x[[1]]
  for(i in 2:length(x)-1){
    xi <- as.vector(outer(xi,x[[i+1]],paste))
  }
  lapply(strsplit(xi,' '),function(x){sort(as.numeric(x))})
}
make.constraints_byigraph <- function(mat,score){
  circles <- get_all_circles(mat,score)
  paths <- lapply(circles,function(circlei){
    lapply(lapply(2:length(circlei)-1,function(i){circlei[i,i+1]}),function(pathi){
      which(score[,ncol(mat)+1]==pathi[1] & score[,pathi[2]]==1)
    })
  })
  do.call(rbind,lapply(unique(do.call(c,lapply(paths,pastes))),function(x){
    tmp <- rep(0,nrow(score))
    tmp[x] <- 1
    tmp
  }))
}

graph_to_data_frame <- function(g1){
  gdf <- lapply(1:length(V(g1)),function(i){
    to=names(which(g1[i]>0))
    if(length(to)>0){
      return(cbind(from=names(V(g1))[i],to=to))
    } else {
      return(NULL)
    }
  })
  do.call(rbind,gdf)
}
align_network <- function(net1,net2,cut2net1=F){
  
  net1 <- rbind(net1,matrix(0,nrow=ncol(net1) - nrow(net1),ncol=ncol(net1)))
  rownames(net1) <- colnames(net1)
  
  if(cut2net1){
    net2 <- net2[rownames(net2)%in%rownames(net1),colnames(net2)%in%colnames(net1)]
  }
  # if(is.character(all.equal(dimnames(net1),dimnames(net2)))){
    # print('###Warning, dimnames not equal###')
  # }
  g1 <- graph_from_adjacency_matrix(t(net1)); g2 <- graph_from_adjacency_matrix(t(net2))
  if(sum(net1)==0){
    g1.match=0
  } else {
    g1.match <- apply(graph_to_data_frame(g1),1,function(x){
      nfrom=V(g2)[names(V(g2))==x[1]]
      nto=V(g2)[names(V(g2))==x[2]]
      length(all_shortest_paths(g2,nfrom,nto)$res)
      # all_shortest_paths(g2,nfrom,nto)
    })  
  }
  # net2 <- net2[rownames(net2)%in%rownames(net1),colnames(net2)%in%colnames(net1)]
  c(aligned=sum(g1.match>0),net1=sum(net1),net2=sum(net2),pd=sum(g1.match>0)/sum(net1),tp=sum(g1.match>0)/sum(net2))
}

#############################
# Process data
#############################

setwd('C:\\Users\\zhu2\\Documents\\getpathway')
load("~/getpathway/gene39761.rda.RData")
load("C:/Users/zhu2/Documents/getpathway/ptwmap2.rda")

#Process data

expinpath <- expinpath[-3]
inputs <- lapply(unique(pathlist[,1]),function(grpi){
  paths <- pathlist[pathlist[,1]==grpi,2]
  input <- expinpath[names(expinpath)%in%paths]
  print(list(groupname=grpi,paths=names(input)))
  input
})
names(inputs) <- unique(pathlist[,1])
inputs <- inputs[sapply(inputs,length)>0]

pcas <- lapply(inputs,function(x){
  x <- do.call(cbind,x)
  x <- qpca(x,min(30,which(pca(x)$prop>=.9)[1]))
  x$X[,1:which(x$prop>=0.9)[1],drop=F]
})

#Load network
setwd("C:/Users/zhu2/Documents/getpathway/model20170215")
load(file='rlt_p2pinp.rda')
load(file='rlt_p2pcp.rda')
load(file='pathwaymap.rda')

#############################
# Modeling
#############################

processx <- function(x,y){
  temp <- matrix(0,nrow=nrow(y)*ncol(y),ncol=ncol(x)*ncol(y))
  for(i2 in 1:ncol(y)){
    temp[(i2-1)*nrow(y)+1:nrow(y),(i2-1)*ncol(x)+1:ncol(x)] <- scale(x)
  }
  temp
}

equationj2 <- function(j,i.y,lambda=0.5){
  i.x <- which(rlt_p2pcp[[2]][i.y,,drop=F]==1)
  y.i <- lapply(rlt_p2pinp[[1]][[i.y]],function(x){x})
  x.i <- lapply(do.call(c,rlt_p2pinp[[1]][i.x]),function(x){x})
  names(x.i) <- unlist(lapply(inputs[i.x],names))
  
  y <- y.i[[j]]
  x <- c(y.i[-j],x.i)
  #process data
  x.p <- lapply(x,function(xi){
    processx(xi,y)
  })
  for(i2 in 1:length(x.p)){colnames(x.p[[i2]]) <- paste(names(x.p)[i2],1:ncol(x.p[[i2]]),sep='_')}
  y <- cbind(scale(as.vector(y)))
  #group lasso model
  gc()
  lambda <- lambdamax(x=cbind(1,do.call(cbind,x.p)),y=y, 
                      index=c(NA,rep(1:length(x),sapply(x.p,ncol))), 
                      penscale = sqrt, model = LinReg(),
                      center=TRUE,standardized=TRUE) * 0.5^(1/lambda-1)
  fit <- grplasso(x=cbind(1,do.call(cbind,x.p)),y=y,
                  index=c(NA,rep(1:length(x),sapply(x.p,ncol))),lambda=lambda,model=LinReg(),
                  penscale = sqrt,
                  control = grpl.control(update.hess = "lambda", trace = 0))
  temp <- rep(0,length(c(y.i,x.i)))
  temp[-j] <- tapply(coef(fit)[-1]!=0,rep(1:length(x),sapply(x.p,ncol)),function(x){any(x!=0)})
  names(temp) <- c(names(y.i),names(x.i))
  temp
}

model_connectout <- function(i.y,lambda=0.5){
  j.len <- length(rlt_p2pinp[[1]][[i.y]])
  modeli <- lapply(1:j.len,function(j){
    # print(j)
    equationj2(j,i.y=i.y,lambda=lambda)
  })
  temp <- do.call(rbind,modeli)
  dimnames(temp) <- list(
    names(modeli[[1]])[1:nrow(temp)],names(modeli[[1]])
  )
  temp[1:nrow(temp),1:nrow(temp)] <- rlt_p2pinp[[2]][[i.y]]
  temp
}

# modeling

load('temp.rda')
# test <- lapply(1:length(rlt_p2pinp[[1]]),function(i){
#   print(i)
#   gc()
#   temp <- lapply((6:14)/20,function(lambda){
#     print(lambda)
#     model_connectout(i,lambda)
#   })
#   # temp <- model_connectout(i,lambda=0.6)
#   temp
# })

# test2 <- lapply(1:length(rlt_p2pinp[[1]]),function(i){
#   print(i)
#   gc()
#   temp <- lapply((8:10)/10,function(lambda){
#     print(lambda)
#     model_connectout(i,lambda)
#   })
#   # temp <- model_connectout(i,lambda=0.6)
#   temp
# })
# test3 <- lapply(1:44,function(i){
#   c(test[[i]],test2[[i]])
# })

#############################

# i <- 0
test_vali <- lapply(test3,function(x){
  # print(i <<- i+1)
  # j <<- 0
  t(sapply(x,function(xi){
    # print(j <<- j+1)
    round(align_network(xi,pthnet),4)
  }))
})

i <- 0
print(i <<- i+1); test_vali[[i]]; sum(test[[i]][[1]][,1:nrow(test[[i]][[1]])]);plotnet(test[[i]][[1]][,1:nrow(test[[i]][[1]])])

#############################

ind <- c(12,1,1,4,6,1,8,7,2,3,12,12,12,12,11,12,12,10,7,11,12,5,4,12,12,10,12,12,9,10,11,12,12,12,11,12,11,12,12,12,12,12,12,12)
test <- lapply(1:44,function(i){test3[[i]][[ind[i]]]})
test2 <- matrix(0,length(unlist(lapply(rlt_p2pinp[[1]],names))),length(unlist(lapply(rlt_p2pinp[[1]],names))),
                dimnames=list(unlist(lapply(rlt_p2pinp[[1]],names)),unlist(lapply(rlt_p2pinp[[1]],names))))
for(i in 1:length(test)){
  test2[match(rownames(test[[i]]),rownames(test2)),match(colnames(test[[i]]),colnames(test2))] <- test[[i]]
}
rlt_p2pfinal <- list(test,test2)
save(rlt_p2pfinal,file='rlt_p2pfinal.rda')

##############################
# Validation
##############################

align_network(test2,pthnet)
