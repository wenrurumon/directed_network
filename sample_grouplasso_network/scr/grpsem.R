
##############################
# Macro
##############################

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

pca2 <- function(Y,prop=0.8){
  Y.rank <- which(qpca(Y,rank=0,ifscale=TRUE)$prop>=sqrt(prop))[1]
  rlt <- qpca(Y,rank=Y.rank,ifscale=TRUE)
  rlt <- rlt$X[,1:which(rlt$prop>=sqrt(prop))[1],drop=F]
  rlt
}

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
grpsem <- function(Y,lambda=.5,outaslist=FALSE){
  models <- lapply(1:length(Y),function(i){
  	print(i)
    data <- data4sem(y=Y[[i]],x=Y[-i])
    coef <- rep(0,length(Y))
    coef[-i] <- quick_grplasso(y=data$y,x=cbind(1,data$x),index=c(NA,data$index),lambda)
    as.numeric(coef)
  })
  if(!outaslist){
    return(do.call(rbind,models))
  } else {
    return(models)
  }
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
                      print(i)
                      cbind(subpath(x.graph[i,],max.parent=max.parent),i)
                    }))
  x.cost <- sapply(1:nrow(x.path),function(i){
    print(i)
    mape(Y[[x.path[i,ncol(x.path)]]],
         Y[which(x.path[i,-ncol(x.path)]>0)])
  })
  return(list(x.graph=x.graph,x.path=x.path,x.cost=x.cost))
}
#
model <- function(Y,x.graph=NULL,prop=NULL,lambda=.5,max.parent=3){
  gc()
  if(length(Y)==1){
    rlt2 <- matrix(0,1,1,dimnames=list(names(Y),names(Y)))
  } else {
    if(!is.null(prop)){
      Y <- lapply(Y,pca2,prop=prop)
    }
    data1 <- ip_data(Y,x.graph,lambda,max.parent)
    rlt1 <- ip_zy(data1$x.graph,data1$x.path,data1$x.cost)
    data2 <- ip_data(Y,rlt1,lambda,max.parent)
    rlt2 <- ip_hzx(data2$x.graph,data2$x.path,data2$x.cost)
    dimnames(rlt2) <- list(names(Y),names(Y))  
  }
  return(rlt2)
}
