
rm(list=ls())

##############################
# Result Validation
##############################

library(igraph)
library(WGCNA)
library(grplasso)
setwd('C:\\Users\\zhu2\\Documents\\getpathway\\model20170126')
# for(i in dir(pattern='rda')){load(i)}

load('pathwaymap.rda')
summary_signalnet <- function(mat){
  g <- graph_from_adjacency_matrix(t(mat))
  temp <- lapply(1:length(pathwaymap),function(i){
    pathi <- unique(pathwaymap[[i]])[!grepl(names(pathwaymap)[i],unique(pathwaymap[[i]]))]
    node.from <- V(g)[grep(names(pathwaymap)[i],names(V(g)))]
    node.to <- V(g)[grep(paste(pathwaymap[[i]],collapse="|"),names(V(g)))]
    paths <- all_shortest_paths(g,node.from,node.to)$res
    list(paths=paths,ends=unique(sapply(paths,function(x){x[length(x)]})))
  })
  list(sum(sapply(temp,function(x){length(x$ends)})),
       sum(sapply(pathwaymap,length)),
       sum(sapply(temp,function(x){length(x$ends)}))/sum(sapply(pathwaymap,length)))
}
summary_signalnet_local <- function(mat){
  if(!is.matrix(mat)){
    return(rep(NA,5))
  }
  g <- graph_from_adjacency_matrix(t(mat))
  pathwaymap <- pathwaymap[names(pathwaymap)%in%colnames(mat)]
  pathwaymap <- lapply(pathwaymap,function(x){x[x%in%colnames(mat)]})
  pathwaymap <- pathwaymap[sapply(pathwaymap,length)>0]
  if(length(pathwaymap)==0){
    return(c(NA,NA,sum(mat),NA,0))
  }
  
  temp <- lapply(1:length(pathwaymap),function(i){
    pathi <- unique(pathwaymap[[i]])[!grepl(names(pathwaymap)[i],unique(pathwaymap[[i]]))]
    node.from <- V(g)[grep(names(pathwaymap)[i],names(V(g)))]
    node.to <- V(g)[grep(paste(pathwaymap[[i]],collapse="|"),names(V(g)))]
    paths <- all_shortest_paths(g,node.from,node.to)$res
    list(paths=paths,ends=unique(sapply(paths,function(x){x[length(x)]})))
  })
  
  pathalign=sum(sapply(temp,function(x){length(x$ends)}))
  pathkegg=sum(sapply(pathwaymap,length))
  pathfound=sum(mat)
  c(pathalign,pathkegg,pathfound,pathalign/pathkegg,pathalign/pathfound)

}

summary_signalnet2 <- function(mat){
  g <- graph_from_adjacency_matrix(t(mat))
  mat3 <- mat2 <- matrix(0,nrow(mat),ncol(mat),dimnames=dimnames(mat))
  temp <- list()
  for(i in 1:length(pathwaymap)){
    pathi <- unique(pathwaymap[[i]])[!grepl(names(pathwaymap)[i],unique(pathwaymap[[i]]))]
    if(length(pathi)==0){next}
    node.from <- V(g)[grep(names(pathwaymap)[i],names(V(g)))]
    node.to <- V(g)[grep(paste(pathwaymap[[i]],collapse="|"),names(V(g)))]
    paths_all <- all_simple_paths(g,node.from,node.to)
    for(pathj in paths_all){
      for(j in 1:(length(pathj)-1)){
        mat2[grep(names(pathj[j+1]),colnames(mat2)),grep(names(pathj[j]),colnames(mat2))] <- 1
      }
    }
    paths_short <- all_shortest_paths(g,node.from,node.to)$res
    paths_short <- paths_short[sapply(paths_short,length)>1]
    for(pathj in paths_short){
      for(j in 1:(length(pathj)-1)){
        mat3[grep(names(pathj[j+1]),colnames(mat3)),grep(names(pathj[j]),colnames(mat3))] <- 1
      }
    }    
    temp[[i]] <- list(paths_all=paths_all,paths_short=paths_short,ends=unique(sapply(paths_short,function(x){x[length(x)]})))
  }
  dp <- sum(sapply(temp,function(x){length(x$ends)}))/sum(sapply(pathwaymap,length))
  fdr_all <- (sum(mat)-sum(mat2))/sum(mat)
  fdr_short <- (sum(mat)-sum(mat3))/sum(mat)
  list(path_summary=temp,g.all=mat2,g.short=mat3,kpi=c(dp=dp,fdr_all=fdr_all,fdr_short=fdr_short))
}

cormat <- function(x,y){
  cor(as.vector(x),as.vector(y))
}
cormatlist <- function(X){
  names(X) <- 1:length(X)
  X <- X[sapply(X,is.matrix)]
  sapply(X,function(x){
    sapply(X,function(y){
      cormat(x,y)
    })
  })
}

##############################
# Macros_processing
##############################

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
plotnet <- function(x,mode='directed',main=NULL){
  diag(x) <- 0
  plot(graph_from_adjacency_matrix(t(x),mode=mode),
       edge.arrow.size=.3,
       vertex.size=.1,
       vertex.label.cex=1.2,
       edge.width=.1,
       main=main)
}

##############################
# Macro_grplasso_sem
##############################

library(grplasso)
#expand variable x to equation set
expandx <- function(x,n){
  temp <- lapply(1:n,function(i){
    temp <- matrix(0,nrow=nrow(x)*n,ncol=ncol(x))
    temp[1:nrow(x)+(i-1)*nrow(x),] <- as.matrix(x)
    temp
  })
  do.call(cbind,temp)
}
#prepare data for equation set
data4sem <- function(y,x){
  x2 <- lapply(x,function(xi){
    expandx(xi,ncol(y))
  })
  list(
    y=as.vector(y),
    x=do.call(cbind,x2),
    index=rep(1:length(x2),sapply(x2,ncol))
  )
}
#quick group lasso sturcture
quick_grplasso <- function(y,x,index,lambda){
  lambda <- lambdamax(x=as.matrix(x),y=as.matrix(y),index=index, 
                      penscale = sqrt, model = LinReg(),
                      center=TRUE,standardized=TRUE) * 0.5^(1/lambda-1)
  fit <- grplasso(x=as.matrix(x),y=as.matrix(y),index=index,lambda=lambda,model=LinReg(),
                  penscale = sqrt,
                  control = grpl.control(update.hess = "lambda", trace = 0))
  tapply(coef(fit),index,function(x){any(x!=0)})
}
#raw group lasso sem within Ys
grpsem <- function(Y,lambda=.5){
  models <- lapply(1:length(Y),function(i){
    data <- data4sem(y=Y[[i]],x=Y[-i])
    coef <- rep(0,length(Y))
    coef[-i] <- quick_grplasso(y=data$y,x=cbind(1,data$x),index=c(NA,data$index),lambda)
    as.numeric(coef)
  })
  do.call(rbind,models)
}
#generate subpaths for pathi with specific max.parent setup
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
#mape for the equation
mape <- function(y,x){
  if(length(x)==0){
    return(mean((scale(y))^2))
  } else {
    x <- do.call(cbind,x)
    return(mean((lm(scale(y)~scale(x))$residual)^2)  )
  }
}
#Integer Programming
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
#
ip_data <- function(Y,x.graph=NULL,lambda=.5,max.parent=3){
  if(is.null(x.graph)){
    x.graph <- grpsem(Y,lambda=lambda)
  }
  x.path <- do.call(rbind,
                    lapply(1:nrow(x.graph),function(i){
                      cbind(subpath(x.graph[i,],max.parent=max.parent),i)
                    }))
  x.cost <- sapply(1:nrow(x.path),function(i){
    mape(Y[[x.path[i,ncol(x.path)]]],
         Y[which(x.path[i,-ncol(x.path)]>0)])
  })
  return(list(x.graph=x.graph,x.path=x.path,x.cost=x.cost))
}
#
model <- function(Y,x.graph=NULL,lambda=.5,max.parent=3){
  gc()
  if(length(Y)==1){
    rlt2 <- matrix(0,1,1,dimnames=list(names(Y),names(Y)))
  } else {
    data1 <- ip_data(Y,x.graph,lambda,max.parent)
    rlt1 <- ip_zy(data1$x.graph,data1$x.path,data1$x.cost)
    data2 <- ip_data(Y,rlt1,lambda,max.parent)
    rlt2 <- ip_hzx(data2$x.graph,data2$x.path,data2$x.cost)
    dimnames(rlt2) <- list(names(Y),names(Y))  
  }
  return(rlt2)
}
test <- function(Y){
  t0 <- try(model(Y,lambda=0.5,max.parent=3))
  t01 <- try(model(Y,lambda=0.5,max.parent=2))
  t1 <- try(model(Y,lambda=0.55,max.parent=3))
  t2 <- try(model(Y,lambda=0.55,max.parent=2))
  t3 <- try(model(Y,lambda=0.6,max.parent=3))
  t4 <- try(model(Y,lambda=0.6,max.parent=2))
  t5 <- try(model(Y,lambda=0.65,max.parent=3))
  t6 <- try(model(Y,lambda=0.65,max.parent=2))
  t7 <- try(model(Y,lambda=0.7,max.parent=3))
  t8 <- try(model(Y,lambda=0.7,max.parent=2))
  t9 <- try(model(Y,lambda=0.75,max.parent=3))
  t10 <- try(model(Y,lambda=0.75,max.parent=2))
  list(t0,t01,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10)
}

##############################
# Process data
##############################

setwd('C:\\Users\\zhu2\\Documents\\getpathway')
load("~/getpathway/gene39761.rda.RData")
load("C:/Users/zhu2/Documents/getpathway/ptwmap2.rda")

pcas <- lapply(expinpath,function(x){
  if(ncol(x)==1){
    return(x)
  }else{
    x <- qpca(x,which(pca(x)$prop>=0.9)[1])
    scale(x$X[,1:which(x$prop>=0.9)[1],drop=F])
  }  
})
for(i in 1:length(pcas)){
  colnames(pcas[[i]]) <- paste(ptwmap[match(names(pcas),ptwmap[,2]),1][i],1:ncol(pcas[[i]]),sep="_")
}
inputs <- lapply(unique(pathlist[,1]),function(grpi){
  paths <- pathlist[pathlist[,1]==grpi,2]
  input <- pcas[names(pcas)%in%paths]
  print(list(i=i,groupname=grpi,paths=names(input)))
  input
})
names(inputs) <- unique(pathlist[,1])
inputs <- inputs[sapply(inputs,length)>0]

#model and validation
load('model20170126/test_rlt.rda')
# test_rlt <- lapply(1:length(inputs),function(i){
#   print(i)
#   test(inputs[[i]])
# })

test_vali <- lapply(test_rlt,function(x){
  sapply(x,summary_signalnet_local)
})

i <- 0

print(i<-i+1); t(test_vali[[i]])
par(mfrow=c(1,1))
test_rlti <- test_rlt[[i]]
corrplot::corrplot(cormatlist(test_rlti))

print(i)
t(test_vali[[i]])

par(mfrow=c(3,4))
for(j in 1:12){
  try(plotnet(test_rlti[[j]],main=j))
}

ind <- c(11,11,1,5,1,9,7,9,9,1,1,1,1,7,7,7,1,5,5,7,1,11,8,5,4,1,4,9,7,11,3,1,1,9,10,11,7,10,11,11,12,11,11,9)
rlt <- lapply(1:44,function(i){test_rlt[[i]][[ind[[i]]]]})
rlt_p2pinp <- list(inputs=inputs,rlt=rlt)

par(mfrow=c(2,2))
for(i in 1:44){
  plotnet(rlt[[i]],main=names(inputs)[i])
}
