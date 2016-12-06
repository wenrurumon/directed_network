
rm(list=ls())

#################################
# SEM-Preprocess
#################################

#bind to array
abind <- function(...){
  x <- list(...)
  rlt <- array(NA,dim=c(dim(x[[1]]),length(x)))
  for (i in 1:length(x)){rlt[,,i]<-x[[i]]}
  return(rlt)
}

#ADJ aggregation
adj.group <- function(dag,Y.group){
  dag.group <- matrix(0,max(Y.group),max(Y.group))
  for(i in 1:nrow(dag.group)){
    for(j in 1:ncol(dag.group)){
      dag.group[j,i] <- sum(apply(dag[Y.group==i,Y.group==j,drop=F],1,max))
    }
  }
  return(dag.group)
}

#single directed score
mat.sds <- function(dag.group){
  dag_sel <- apply(abind(dag.group,t(dag.group)),1:2,function(x){which(x==max(x))[1]})
  dag.group[dag_sel==2] <- 0; diag(dag.group) <- 0
  dag.group
}

#################################
# SEM-Undirected Structure
#################################

#SEM L1
sem_l1 <- function(Y,lambda=0.1,times=10){
  adjs <- lapply(1:times,function(i){
    Y <- Y[sample(1:nrow(Y),nrow(Y)*2/3),]
    adj <- do.call(rbind,lapply(1:ncol(Y),function(i){
      slimi <- slim(X=Y[,-i],Y=Y[,i],lambda=lambda,
                    rho=1,method='lasso',verbose=FALSE)
      temp <- rep(FALSE,ncol(Y))
      temp[-i][which(slimi$beta!=0)] <- TRUE
      temp
    }))
  })
  adj <- apply(do.call(abind,adjs),1:2,mean)
  return(adj)
}

#SEM Group Lasso
sem_grplasso <- function(Y
                         ,Y.group=rep(1:ncol(Y))
                         ,Y.prop=(1/tapply(Y.group,Y.group,length))[Y.group]
                         ,lambda=.1,times=10){
  adjs <- lapply(1:times,function(i){
    Y <- Y[sample(1:nrow(Y),nrow(Y)*2/3),]
    adj <- do.call(rbind,lapply(1:ncol(Y),function(i){
      lambda <- lambdamax(x=cbind(1,Y[,-i]),y=Y[,i], 
                          index=c(NA,Y.group[-i]), 
                          penscale = sqrt, model = LinReg(),
                          center=TRUE,standardized=TRUE) * 0.5^(1/lambda-1)
      fit <- grplasso(x=cbind(1,Y[,-i]),y=Y[,i],
                      index=c(NA,Y.group[-i]),lambda=lambda,model=LinReg(),
                      penscale = sqrt,
                      control = grpl.control(update.hess = "lambda", trace = 0))
      temp <- rep(0,ncol(Y))
      temp[-i] <- coef(fit)[-1]
      temp!=0
    }))
  })
  adj <- apply(do.call(abind,adjs),1:2,mean)
  return(adj)
}

#SEM L1 with parameter estimation and X connected
sem_l1_YX <- function(Y,X=NULL,adj=NULL,lambda=0.1,times=10,stability=0.8){
  #Lasso Network cross Y
  if(is.null(adj)|is.null(X)){
    adjs <- lapply(1:times,function(i){
      Y <- Y[sample(1:nrow(Y),nrow(Y)*2/3),]
      adj <- do.call(rbind,lapply(1:ncol(Y),function(i){
        slimi <- slim(X=Y[,-i],Y=Y[,i],lambda=lambda,
                      rho=1,method='lasso',verbose=FALSE)
        temp <- rep(FALSE,ncol(Y))
        temp[-i][which(slimi$beta!=0)] <- TRUE
        temp
      }))
    })
    adj <- apply(do.call(abind,adjs),1:2,mean)
    dimnames(adj) <- list(colnames(Y),colnames(Y))
  }
  
  #Parameter Estimation
  if(is.null(X)){
    model <- lapply(1:nrow(adj),function(i){
      yi <- Y[,adj[i,]>=stability,drop=F]
      if(length(yi)==0){return(NULL)}
      yi.coef <- ginv(t(yi)%*%yi) %*% t(yi) %*% Y[,i]  
      yi.sigma <- (1/nrow(yi))*sum((Y[,i]-yi%*%yi.coef)^2)
      yi.SIGMA <- yi.sigma * ginv(t(yi)%*%yi)
      yi.pvalue <- pchisq(as.vector(yi.coef^2)/diag(yi.SIGMA),df=1,lower.tail=FALSE)
      data.frame(y=colnames(Y)[i],x=colnames(yi),coef=as.vector(yi.coef),pvalue=yi.pvalue)
    })
    model <- do.call(rbind,model)
    #Model Output
    return(list(adj=adj,model=model))
  } 
  
  #Connection with X
  adjs <- lapply(1:times,function(timei){
    if(times>1){
      sampleset <- sample(1:nrow(Y),nrow(Y)/2)
      Y <- Y[sampleset,,drop=F]; X <- X[sampleset,,drop=F]
    }
    Xinv<-ginv(t(X)%*%X)
    Yhat<-X%*%Xinv%*%t(X)%*%Y
    adj <- t(sapply(1:ncol(Y),function(k){
      Xk <- as.matrix(cbind(Yhat[,-k],X))
      outl1=slim(X=Xk[,apply(Xk,2,var)!=0],Y=Y[,k],lambda=lambda,
                 rho=1,method = "lasso",verbose=FALSE)
      temp <- rep(FALSE,ncol(Y)+ncol(X))
      temp[-k][apply(Xk,2,var)!=0] <- outl1$beta!=0
      temp
    }))
    return(adj)
  })
  adjYX <- apply(do.call(abind,adjs),1:2,mean)
  dimnames(adjYX) <- list(colnames(Y),c(colnames(Y),colnames(X)))
  # adjYX[,1:ncol(Y)] <- apply(abind(adj,adjYX[,1:ncol(Y)]),1:2,max)
  adjYX[,1:ncol(Y)] <- adj
  
  #Parameter Estimatin
  model <- do.call(rbind,lapply(1:nrow(adjYX),function(i){
    yi <- cbind(Y,X)[,adjYX[i,]>=stability,drop=F]
    if(length(yi)==0){return(NULL)}
    yi.coef <- ginv(t(yi)%*%yi) %*% t(yi) %*% Y[,i]  
    yi.sigma <- (1/nrow(yi))*sum((Y[,i]-yi%*%yi.coef)^2)
    yi.SIGMA <- yi.sigma * ginv(t(yi)%*%yi)
    yi.pvalue <- pchisq(as.vector(yi.coef^2)/diag(yi.SIGMA),df=1,lower.tail=FALSE)
    data.frame(y=colnames(Y)[i],x=colnames(yi),coef=as.vector(yi.coef),pvalue=yi.pvalue)
  }))
  
  #Model Output
  return(list(adj=adjYX,model=model))
}

#grouplasso model with L1 regulazation applied
sem_grplasso2 <- function(Y
                          ,Y.group=rep(1:ncol(Y))
                          ,Y.prop=(1/tapply(Y.group,Y.group,length))[Y.group]
                          ,lambda1=.5,lambda2=.3,times=10,stability=.8){
  #Group Lasso Network
  sem1 <- sem_grplasso(Y,Y.group,lambda=lambda1,times=times)
  adj <- (sem1>=stability)
  #L1 Penalty
  adjs <- lapply(1:times,function(i){
    Y <- Y[sample(1:nrow(Y),nrow(Y)*2/3),]
    adj <- do.call(rbind,lapply(1:nrow(adj),function(j){
      temp <- adj[j,]
      if(sum(temp)==0){return(temp)}
      Yj <- Y[,j,drop=F]
      Xj <- Y[,temp,drop=F]
      slimi <- slim(X=Xj,Y=Yj,lambda=lambda2,rho=1,method='lasso',verbose=FALSE)
      temp[temp] <- (slimi$beta!=0)
      return(temp)
    }))
    return(adj)
  })
  sem2 <- apply(do.call(abind,adjs),1:2,mean)
  #Result
  list(sem_l1=sem1,sem_grplasso=sem2)
}


#####################################################
# Preprocess
#####################################################

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
       edge.arrow.size=.1,
       vertex.size=3,
       vertex.label.cex=1,
       edge.width=.1)
}

###########################################################
# Input data
###########################################################

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
library(flare)
library(grplasso)
library(data.table)
library(dplyr)
library(igraph)

###########################################################
# Test
###########################################################

Y.exp <- expinpath[-3]
pathway.group <- tapply(pathlist[,2]%in%names(Y.exp),pathlist[,1],sum)
pathway.group <- pathway.group[pathway.group>0]
pathway.map <- lapply(names(pathway.group),function(grpi){
  paths <- pathlist[pathlist[,1]==grpi,2]
  names(Y.exp )[names(Y.exp )%in%paths]
})
Y.exp <- Y.exp[names(Y.exp)%in%unlist(pathway.map)]

test <- function(i){
  #Debug
  print(i)
  lambda = 0.3
  stability = 0.8
  #
  Y <- Y.exp[[i]]
  system.time(Y.sem <- sparse_2sem(Y,lambda=lambda,times=1))
  system.time(Y.cnif <- try(CNIF(data=Y,init.adj=(Y.sem[[1]]>=stability),max_parent=3)))
  if(!is.list(Y.cnif)){
    system.time(Y.cnif <- try(CNIF(data=Y,init.adj=(Y.sem[[1]]>=stability),max_parent=2)))
  }
  system.time(Y.rlt <- sparse_2sem(Y,lambda=lambda,Y.fixed=Y.cnif))
  plotnet(Y.rlt[[1]],'directed')
  list(data=Y,adj=Y.rlt)
}
# temp <- test(115)
# plotnet(temp$adj[[1]])

# load('expr_network_20161206.rda')
# length(rlt)

rlt <- list()
for(i in (length(rlt)+1):length(Y.exp)){
  temp <- try(test(i))
  rlt[[i]] <- temp
  save(rlt,file='expr_network_20161206.rda')
}

#Validation
load('expr_network_20161205.rda')
