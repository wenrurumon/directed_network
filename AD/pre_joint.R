
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
setwd('e:\\uthealth\\signaling\\codes\\')
source('sparse_2sem_final.R')
source('local_cnif_macro.R')
source('CNIF.R')
sourceCpp("score_function_regression.cpp")
sourceCpp("simple_cycle.cpp")
sourceCpp("initial_sem.cpp")
setwd('e:\\uthealth\\getpathway')
load("e:\\uthealth\\getpathway/gene39761.rda.RData")
load('e:\\uthealth\\getpathway/model20170215/rlt_p2pinp.rda')
library(flare)
library(grplasso)
library(data.table)
library(dplyr)
library(igraph)

###########################################################
# Get genelist and expression data
###########################################################

load('e:\\uthealth\\getpathway\\AD_data.rda')

# library(rvest)
# 
# kegg <- read_html('http://www.kegg.jp/dbget-bin/www_bget?hsa05010')
# kegg <- kegg %>% html_nodes('tr:nth-child(9) .td30 div') %>% html_text()
# genelist1 <- substr(kegg[1:171*2],1,regexpr(';',kegg[1:171*2])-1)
# 
# kegg <- readLines('http://www.kegg.jp/kegg-bin/show_pathway?map=hsa05010&show_description=show')
# kegg <- kegg[(grep('mapdata',kegg)[2]+1):(grep('</map>',kegg)-1)]
# kegg <- kegg[!grepl('onmouseover',kegg)]
# kegg <- substr(kegg,regexpr('title',kegg)+7,nchar(kegg))
# kegg <- substr(kegg,1,regexpr('\"',kegg)-1)
# kegg <- strsplit(gsub(',','',kegg),' ')
# genelist2 <- lapply(kegg,function(x){
#   x <- x[1:(length(x)/2)*2]
#   substr(x,2,nchar(x)-1)
# })
# genelist2 <- unique(genelist2)
# 
# expinpath <- do.call(cbind,lapply(expinpath,scale))
# exp1 <- expinpath[,match(genelist1,colnames(expinpath))]
# exp2 <- lapply(genelist2,function(x){
#   x <- expinpath[,match(x,colnames(expinpath)),drop=F]
#   x[is.na(x)] <- 0
#   x[,apply(x,2,var)>0,drop=F]
# })
# exp2 <- exp2[sapply(exp2,ncol)>0]
# exp2.gene <- sapply(exp2,function(x){colnames(x)[1]})
# exp2.prop <- sapply(exp2,function(x){qpca(x)$prop[1]})
# exp2 <- sapply(exp2,function(x) qpca(x)$X[,1,drop=F])
# colnames(exp2) <- exp2.gene
# save(genelist1,genelist2,exp1,exp2,file='AD_data.rda')

###########################################################
# expression network
###########################################################

#For exp2

# Y.sem <- sparse_2sem(exp2,lambda=0.08,times=10)
Y.sem <- sparse_2sem(exp2,lambda=0.1,times=10)
plotnet(Y.sem[[1]]>=.8)
# Y.cnif <- CNIF(data=exp2,init.adj=Y.sem[[1]]>=.8,max_parent=2)
Y.cnif <- CNIF(data=exp2,init.adj=Y.sem[[1]]>=.8,max_parent=3)
plotnet(Y.cnif,'directed')
Y.sem <- sparse_2sem(exp2,Y.fixed=Y.cnif)

#For Adni

load('E:\\uthealth\\joint_model\\adni\\adni2.rda')
adni2 <- adni2[,colnames(adni2)%in%colnames(exp2),drop=F]
# adni.sem <- sparse_2sem(adni2,lambda=0.08,times=10)
adni.sem <- sparse_2sem(adni2,lambda=0.1,times=10)
plotnet(adni.sem[[1]]>=0.8)
# adni.cnif <- CNIF(data=adni2,init.adj=adni.sem[[1]]>=.8,max_parent=2)
adni.cnif <- CNIF(data=adni2,init.adj=adni.sem[[1]]>=.8,max_parent=3)
plotnet(adni.cnif,'directed')

# setwd('E:\\uthealth\\20171216\\expression network')
# write.csv(Y.sem$eq_scorenet,'ADrush.csv',row.names=F)
# write.csv(adni.sem$eq_scorenet,'ADadni.csv',row.names=F)

###########################################################
# Overlap
###########################################################

g.rush <- graph_from_adjacency_matrix(t(Y.sem[[1]]),mode='directed')
g.adni <- graph_from_adjacency_matrix(t(adni.sem[[1]]),mode='directed')
m.rush <- Y.sem[[2]]
m.adni <- adni.sem[[2]]
x <- m.adni
x <- lapply(1:nrow(x),function(i){
  x <- x[i,c(2,4)]
  shortest_paths(g.rush,from=x[2],to=x[1])$vpath
})
sum(sapply(lapply(x,unlist),length)>0)/length(x)
x <- m.rush
x <- lapply(1:nrow(x),function(i){
  x <- x[i,c(2,4)]
  shortest_paths(g.adni,from=x[2],to=x[1])$vpath
})
sum(sapply(lapply(x,unlist),length)>0)/length(x)

###########################################################
# Rolling
###########################################################

test <- function(la,ti,mp){
  gc()
  Y.sem <- sparse_2sem(exp2,lambda=la,times=ti)
  plotnet(Y.sem[[1]]>=.8)
  Y.cnif <- CNIF(data=exp2,init.adj=Y.sem[[1]]>=.8,max_parent=mp)
  plotnet(Y.cnif,'directed')
  Y.sem <- sparse_2sem(exp2,Y.fixed=Y.cnif)
  gc()
  adni.sem <- sparse_2sem(adni2,lambda=la,times=ti)
  plotnet(adni.sem[[1]]>=0.8)
  adni.cnif <- CNIF(data=adni2,init.adj=adni.sem[[1]]>=.8,max_parent=mp)
  gc()
  g.rush <- graph_from_adjacency_matrix(t(Y.sem[[1]]),mode='directed')
  g.adni <- graph_from_adjacency_matrix(t(adni.sem[[1]]),mode='directed')
  m.rush <- Y.sem[[2]]
  m.adni <- adni.sem[[2]]
  x <- m.adni
  x <- lapply(1:nrow(x),function(i){
    x <- x[i,c(2,4)]
    shortest_paths(g.rush,from=x[2],to=x[1])$vpath
  })
  print(xadni <- c(length(x),sum(sapply(lapply(x,unlist),length)>0)/length(x)))
  x <- m.rush
  x <- lapply(1:nrow(x),function(i){
    x <- x[i,c(2,4)]
    shortest_paths(g.adni,from=x[2],to=x[1])$vpath
  })
  print(xrush <- c(length(x),sum(sapply(lapply(x,unlist),length)>0)/length(x)))
  return(bk <<- c(adni=xadni,rush=xrush))
}

# test_0.1_10_3 <- test(0.1,10,3)
# adni1   adni2   rush1   rush2 
# 50.0000  0.1800 48.0000  0.1875 
# test(0.1,10,2)
# adni1      adni2      rush1      rush2 
# 48.0000000  0.1458333 47.0000000  0.1489362 
# test(0.08,10,2)
# adni1       adni2       rush1       rush2 
# 91.00000000  0.08791209 58.00000000  0.18965517 
# test(0.13,10,2)
# adni1      adni2      rush1      rush2 
# 17.0000000  0.2941176 26.0000000  0.1923077 

gc();set.seed(1);test(0.11,1,3);bk1 <- bk;bk

3
adni1      adni2      rush1      rush2 
43.0000000  0.1860465 49.0000000  0.1428571 

# 6,0.12,1,13
# adni1      adni2      rush1      rush2 
# 28.0000000  0.1785714 51.0000000  0.1764706 
