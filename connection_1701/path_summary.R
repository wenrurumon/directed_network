
rm(list=ls())
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

g <- graph_from_adjacency_matrix(t(sem_final_signal[[1]]));sum(t(sem_final_signal[[1]]));summary_signalnet(g)
g <- graph_from_adjacency_matrix(t(rlt_p2pfinal[[1]][[18]][,1:24]));sum(rlt_p2pfinal[[1]][[18]][,1:24]);summary_signalnet(g)
g <- graph_from_adjacency_matrix(t(rlt_p2pfinal[[2]]));summary_signalnet(g)

load('pathwaymap.rda')
temp <- t(rlt_p2pfinal[[2]])
g <- graph_from_adjacency_matrix(temp)
temp <- rownames(rlt_p2pfinal[[2]])
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
set.seed(12345); randmat <- lapply(1:10000,function(i){randomDAGmat(24,0.18)})
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

test <- sapply(randmat,function(x){
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

