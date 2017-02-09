
# rm(list=ls())
library(igraph)
library(WGCNA)
library(grplasso)
setwd('C:\\Users\\zhu2\\Documents\\getpathway\\model20170126')
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

# par(mfrow=c(1,2))
# for(i in grep('Signal',names(rlt_p2pfinal[[1]]))){
#   temp <- rlt_p2pfinal[[1]][[i]]
#   plotnet(temp[,1:nrow(temp)])
# }
# 
# g <- graph_from_adjacency_matrix(t(temp[,1:nrow(temp)]))
# lapply(1:5,function(i){all_simple_paths(g,from=V(g)[i])})
# plotnet(temp[,1:nrow(temp)])
# 
# tmp <- matrix(0,0,2)
# for(i in 1:nrow(rlt_p2pfinal[[2]])){
#   for(j in 1:ncol(rlt_p2pfinal[[2]])){
#     if(rlt_p2pfinal[[2]][i,j]>0){
#       tmp <- rbind(tmp,c(rownames(rlt_p2pfinal[[2]])[i],colnames(rlt_p2pfinal[[2]])[j]))
#     }    
#   }
# }
# colnames(tmp) <- c('y','x')
# write.csv(tmp,'p2pfinal.csv',row.names=F)

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
load('pathwaymap.rda')
temp <- t(rlt_p2pfinal[[2]])
g <- graph_from_adjacency_matrix(temp)
temp <- rownames(rlt_p2pfinal[[1]][[18]])
pathwaymap <- pathwaymap[names(pathwaymap)%in%temp]
pathwaymap <- lapply(pathwaymap, function(x){unique(x[x%in%temp])})

summary_signalnet <- function(g){
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

g <- graph_from_adjacency_matrix(t(sem_final_signal[[1]]));sum(t(sem_final_signal[[1]]));summary_signalnet(g)
g <- graph_from_adjacency_matrix(t(rlt_p2pfinal[[1]][[18]][,1:24]));sum(rlt_p2pfinal[[1]][[18]][,1:24]);summary_signalnet(g)
g <- graph_from_adjacency_matrix(t(rlt_p2pfinal[[2]]));summary_signalnet(g)

load('pathwaymap.rda')
temp <- t(rlt_p2pfinal[[2]])
g <- graph_from_adjacency_matrix(temp)
# temp <- rownames(rlt_p2pfinal[[2]])
temp <- rownames(rlt_p2pfinal[[1]][[18]])
pathwaymap <- pathwaymap[names(pathwaymap)%in%temp]
pathwaymap <- lapply(pathwaymap, function(x){unique(x[x%in%temp])})
summary_signalnet(g)

#############################

signal_score <- lapply(rlt_p2pinp[[1]][[18]],function(x){x[[1]]})
datExpr <- do.call(cbind,signal_score)
# net <- blockwiseModules(datExpr,power=6,maxBlockSize = 6000,TOMType='unsigned',minModuleSize = 30,
#                         reassignThreshold = 0,mergeCutHeight = 0.25,
#                         numericLabels = T,pamRespectsDendro = F,saveTOMs = T,saveTOMFileBase = 'AS-green-FPKM-TOM',
#                         verbose=3)
TOM <- TOMsimilarityFromExpr(datExpr,power=6)
dimnames(TOM) <- list(colnames(datExpr),colnames(datExpr))
diag(TOM) <- 0
# cyt <- exportNetworkToCytoscape(TOM,nodeNames=names(datExpr),threshold=0.02)
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

sum(rlt_p2pinp[[2]][[18]][,1:24])
set.seed(12345); randmat <- lapply(1:100,function(i){randomDAGmat(24,0.18)})
randmat <- lapply(randmat,function(x){
  dimnames(x) <- dimnames(rlt_p2pinp[[2]][[18]])
  return(x)
})
randmat <- c(list(rlt_p2pinp[[2]][[18]][,1:24]),randmat)

summary_signalnet <- function(g){
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

# i <- 0
test <- sapply(randmat,function(x){
  # print(i<<-i+1)
  g <- graph_from_adjacency_matrix(t(x))
  rlt <- try(summary_signalnet(g))
  if(!is.list(rlt)) return(NA)
  c(sum(x),unlist(rlt))
})
if(is.list(test)){
  test <- do.call(rbind,test[sapply(test,length)==4])
} else {
  test <- t(test)
}
head(test)
hist(test[-1,4])
test <- cbind(test,test[,2]/test[,1])
head(test)
par(mfrow=c(1,2))
hist(test[-1,4],main = '#aligned paths/#kegg paths'); hist(test[-1,5],main='#aligned paths/#all paths detected')

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
# plotnet(bnnet>=sort(bnnet[bnnet>0],T)[49])
plot(g.bnnet <- graph_from_adjacency_matrix(t(bnnet>=sort(bnnet[bnnet>0],T)[49])))
summary_signalnet(g.bnnet)
sum(bnnet>=sort(bnnet[bnnet>0],T)[49])

#####################################
# Validation
#####################################

vali_om <- summary_signalnet2(rlt_p2pinp[[2]][[18]][,1:24])$kpi
vali_bn <- summary_signalnet2(bnnet>=sort(bnnet[bnnet>0],T)[49])$kpi
vali_rand <- lapply(randmat[-1],function(x){summary_signalnet2(x)$kpi})

##################################

