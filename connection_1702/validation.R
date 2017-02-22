
rm(list=ls())

library(igraph)
library(WGCNA)
library(grplasso)
setwd('C:\\Users\\zhu2\\Documents\\getpathway\\model20170215')
for(i in dir(pattern='rda')){load(i)}

plotnet <- function(x,mode='directed',main=''){
  diag(x) <- 0
  plot(graph_from_adjacency_matrix(t(x),mode=mode),
       edge.arrow.size=.3,
       vertex.size=3,
       vertex.label.cex=1,
       edge.width=.1,
       main=main)
}

plotnet(rlt_p2pfinal[[2]])
par(mfrow=c(1,1))
for(i in grep('Signal',names(rlt_p2pfinal[[1]]))){
  temp <- rlt_p2pfinal[[1]][[i]]
  plotnet(temp[,1:nrow(temp),drop=F],main=names(rlt_p2pfinal[[1]])[i])
}
par(mfrow=c(2,2))
for(i in 1:length(rlt_p2pfinal[[1]])){
  print(i)
  temp <- rlt_p2pfinal[[1]][[i]]
  plotnet(temp[,1:nrow(temp),drop=F],main=names(rlt_p2pfinal[[1]])[[i]])
}

##############################
#Validation Macro

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

align_network2 <- function(net1,net2,cut2net1=F){
  
  net1 <- rbind(net1,matrix(0,nrow=ncol(net1) - nrow(net1),ncol=ncol(net1)))
  rownames(net1) <- colnames(net1)
  
  if(cut2net1){
    net2 <- net2[rownames(net2)%in%rownames(net1),colnames(net2)%in%colnames(net1)]
  }
  g1 <- graph_from_adjacency_matrix(t(net1)); g2 <- graph_from_adjacency_matrix(t(net2))
  if(sum(net1)==0){
    g1.match=NULL
  } else {
    g1.match <- apply(graph_to_data_frame(g1),1,function(x){
      nfrom=V(g2)[names(V(g2))==x[1]]
      nto=V(g2)[names(V(g2))==x[2]]
      all_shortest_paths(g2,nfrom,nto)$res[1]
    })  
  }
  g1.match <- do.call(c,g1.match)
  list(
    g1.match,
    sum(sapply(g1.match,function(x){!is.null(x)})),
    sum(sapply(g1.match,length))
  )
}

round(align_network(net_om <- rlt_p2pinp[[2]][[18]],pthnet,T),4)

#############################

signal_score <- lapply(rlt_p2pinp[[1]][[18]],function(x){x})
datExpr <- do.call(cbind,signal_score)
TOM <- TOMsimilarityFromExpr(datExpr,power=6)
dimnames(TOM) <- list(colnames(datExpr),colnames(datExpr))
diag(TOM) <- 0
cenet <- (TOM>=quantile(TOM,0.95));sum(cenet)
cenet <- apply(apply(cenet,1,function(x){tapply(x,rep(names(signal_score),sapply(signal_score,ncol)),sum)}),1,function(x){tapply(x,rep(names(signal_score),sapply(signal_score,ncol)),sum)})
par(mfrow=c(1,1)); plotnet(cenet>0,'undirected')

##############################

Y <- signal_score
grplasso_model <- function(yi,xi,lambda=.5){
  x=cbind(1,scale(do.call(cbind,xi)))
  y=scale(yi)
  index=c(NA,rep(1:length(xi),sapply(xi,ncol)))
  lambda <- lambdamax(x=x,y=y,index=index, 
                      penscale = sqrt, model = LinReg(),
                      center=TRUE,standardized=TRUE) * 0.5^(1/lambda-1)
  fit <- grplasso(x=x,y=y,index=index,lambda=lambda,model=LinReg(),
                  penscale = sqrt,
                  control = grpl.control(update.hess = "lambda", trace = 0))
  return(fit)
}
equationi <- function(i,lambda=0.8){
  yi <- signal_score[[i]]
  xi <- signal_score[-i]
  fit <- lapply(1:ncol(yi),function(j){
    grplasso_model(yi[,j],xi,lambda=lambda)
  })
  coef.fit <- sapply(fit,function(x){coef(x)[-1]})
  tmp <- rep(0,length(signal_score))
  tmp[-i] <- rowSums(apply(coef.fit,2,function(x){
    tapply(x,rep(names(xi),sapply(xi,ncol)),function(x){any(x!=0)})
  }))>0
  tmp
}
score_grplasso <- do.call(rbind,lapply(1:length(signal_score),equationi))
dimnames(score_grplasso) <- list(names(signal_score),names(signal_score))
plotnet(score_grplasso,'undirected')

#########################
#Random DAG simulation
#########################

# source('http://bioconductor.org/biocLite.R')
# biocLite(c('graph','RBGL','Rgraphviz'))

library(pcalg)
randomDAGmat <- function(n,p){
  showAmat(randomDAG(n,p))==2
}

randmats <- lapply(1:10,function(i){
  set.seed(i); randmat <- lapply(1:1000,function(i){randomDAGmat(24,0.16)})
  randmat <- lapply(randmat,function(x){
    dimnames(x) <- dimnames(rlt_p2pinp[[2]][[18]])
    return(x)
  })
})

i <- 0
simulation <- lapply(randmats,function(randmat){
  print(i <<- i+1)
  t(sapply(randmat,function(x){
    c(align_network(x,pthnet,T),align_network(x,pthnet,F))
  }))
})
simulation2 <- lapply(randmats,function(randmat){
  print(i <<- i+1)
  t(sapply(randmat,function(x){
    unlist(align_network2(x,pthnet,T)[2:3])
  }))
})

simulation3 <- lapply(randmats,function(randmat){
  print(i <<- i+1)
  lapply(randmat,function(x){
    align_network2(x,pthnet,T)[[1]]
  })
})

vali <- function(xi,j=3){
  xi.len <- length(xi)
  xi <- xi[!sapply(xi,is.null)]
  sum(sapply(xi,length)<=j)/xi.len
}
test <- lapply(simulation3,function(x){sapply(x,vali)})
sapply(test,summary)
vali(xi <- align_network2(net_om,pthnet,T)[[1]])

#################################
# DAG with other methodologies
#################################

# sort(cenet,T)[49*2]
# adj <- (cenet>=7)
# adj2 <- do.call(rbind,lapply(1:nrow(adj),function(i){cbind(subpath(adj[i,],max.parent = 2),i)}))
# cost <- mape(adj2,signal_score)
# coe_dag <- ip_zy(adj,adj2,cost)

##############################
#Bayesian Network
library(bnlearn)
# datExpr <- do.call(cbind,lapply(signal_score,function(x) x[,1]))[,1:10]
datExpr <- do.call(cbind,signal_score)
system.time(bn.hc <- hc(as.data.frame(datExpr)))
bnhc3 <- bn.hc[[3]]
bnnet <- matrix(0,nrow=nrow(TOM),ncol=ncol(TOM),dimnames=dimnames(TOM))
for(i in 1:nrow(bnhc3)){
  bnnet[grep(bnhc3[i,2],rownames(bnnet)),grep(bnhc3[i,1],colnames(bnnet))] <- 1
}
bnnet <- apply(apply(bnnet,1,function(x){tapply(x,rep(names(signal_score),sapply(signal_score,ncol)),sum)}),1,function(x){tapply(x,rep(names(signal_score),sapply(signal_score,ncol)),sum)})
diag(bnnet) <- 0
for(i in 1:nrow(bnnet)){
  for(j in 1:ncol(bnnet)){
    if(bnnet[i,j]>bnnet[j,i]){
      bnnet[j,i] <- 0
    } else {
      bnnet[i,j] <- 0
    }
  }
}
bnnet2 <- bnnet>=sort(bnnet[bnnet>0],T)[49]
vali(align_network2(bnnet2,pthnet,T)[[1]])
vali(align_network2(net_om,pthnet,T)[[1]])

#####################################
# Validation
#####################################

test <- sapply(do.call(c,simulation3),vali); summary(test)
vali(align_network2(bnnet2,pthnet,T)[[1]])
vali(align_network2(net_om,pthnet,T)[[1]])

test <- sapply(do.call(c,simulation3),function(x){vali(x,2)}); summary(test)
vali(align_network2(bnnet2,pthnet,T)[[1]],2)
vali(align_network2(net_om,pthnet,T)[[1]],2)

##################################

