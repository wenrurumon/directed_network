
rm(list=ls())
library(igraph)
library(slam)
setwd('C:/Users/zhu2/Documents/getpathway/model20170215/summary')
load('C:/Users/zhu2/Documents/getpathway/model20170215/rlt_p2pinp.rda')

plotg <- function(g,cg=c(.3,1,1,1)){
  v <- names(V(g))
  V(g)$color <- rep('yellow',length(v))
  V(g)$color[grep("d'",v)] <- 'red'
  V(g)$color[grep("ph'",v)] <- 'green'
  V(g)$color[grep("p'",v)] <- 'blue'
  # names(V(g)) <- gsub("d'|p'|ph'|m'|g'","",names(V(g)))
  plot(g,
       edge.arrow.size=cg[1],
       vertex.size=cg[2],
       vertex.label.cex=cg[3],
       edge.width=cg[4])
}
df2mat <- function(gdf){
  v <- unique(as.vector(gdf))
  out <- matrix(0,length(v),length(v),dimnames=list(v,v))
  out2 <- as.matrix(slam::simple_triplet_matrix(match(gdf[,2],v),match(gdf[,1],v),rep(1,nrow(gdf))))
  out[1:nrow(out2),1:ncol(out2)] <- out2
  out
}
mat2df <- function(xmat){
  temp <- matrix(0,0,3)
  for(i in 1:nrow(xmat)){
    for(j in 1:ncol(xmat)){
      if(xmat[i,j]>0){
        temp <- rbind(temp,c(j-1,i-1,xmat[i,j]))
      }
    }
  }
  colnames(temp) <- c('source','target','value')
  data4plot <- list(edge=temp,node=colnames(xmat))
  x.link <- as.data.frame(temp)
  x.node <- colnames(xmat)
  group=apply(sapply(c("ph'","d'","p'","g'","m'"),function(x){grepl(x,x.node)}),1,which)
  group <- c('pheno','disease','pathway','geno','methylation')[group]#group recoding
  x.node <- gsub("ph'|p'|d'|g'|m'","",x.node) #node recoding
  x.node <- data.frame(name=x.node,group=group,size=1)
  temp[,1] <- colnames(xmat)[as.numeric(temp[,1])+1]
  temp[,2] <- colnames(xmat)[as.numeric(temp[,2])+1]
  list(df=temp,link=x.link,node=x.node)
}

############################################################

gdf <- as.matrix(read.csv('forcyt.csv'))
gdf <- rbind(gdf,
             c("g'CREBBP","d'AD"),
             c("m'POU3F2","d'AD"))
gmat <- df2mat(gdf)
pmap <- lapply(rlt_p2pinp[[1]],names)

#select nodes
pathsel <- unlist(pmap[c(18,29,37)])
node.sel <- grepl("ph'|d'|g'|m'",colnames(gmat))|
            grepl(paste(pathsel,collapse='|'),colnames(gmat))
node.sel <- grepl("ph'|d'",colnames(gmat))
node.sel <- grepl("p'",colnames(gmat))
node.sel <- grepl(paste(pathsel,collapse='|'),colnames(gmat))
node.sel <- grepl("ph'|d'",colnames(gmat))|
            grepl(paste(pathsel,collapse='|'),colnames(gmat))|
            grepl('CREBBP|POU',colnames(gmat))
node.sel <- grepl("ph'|d'",colnames(gmat))

#get data4plot
xmat <- gmat[node.sel,node.sel]
# xmat <- gmat
g <- graph_from_adjacency_matrix(t(xmat))
g.cluster <- clusters(g)
xmat <- xmat[which(g.cluster$membership==which(g.cluster$csize==max(g.cluster$csize))),which(g.cluster$membership==which(g.cluster$csize==max(g.cluster$csize)))]
data4plot <- mat2df(xmat)

#############################################################

# Create graph
head(data4plot$link)
head(data4plot$node)
data4plot$node$size <- 
  c(30,5,1,0.01,0.01)[match(data4plot$node$group,c('disease','pheno','pathway','geno','methylation'))]

xlink <- data4plot$link
xnode <- data4plot$node
# xnode$group <- names(pmap)[sapply(paste(xnode$name),function(x){x <- grep(x,pmap);ifelse(length(x)==0,8,x)})]
library(networkD3)
forceNetwork(Links = xlink, Nodes = xnode, Source = "source",
             Target = "target", Value = "value", NodeID = "name",
             Nodesize = "size",linkColour = "#999",
             radiusCalculation = "Math.sqrt(d.nodesize)+6",
             Group = "group", opacity =1,charge=-30, legend = F
             ,zoom=T,opacityNoHove=1)

xlink <- cbind(source=c(0,0,1),target=c(1,2,2))
xlink <- rbind(xlink,xlink+3,xlink+6)
xlink <- data.frame(rbind(xlink,c(3,1),c(7,5),c(9,7),c(10,8),c(9,4),c(10,3)),value=1)
xnode <- data.frame(name=0:10,
      group=rep(c('disease','pheno','pathway','geno','methylation'),c(3,3,3,1,1)),
      size=rep(c(30,5,1,0.01,0.01),c(3,3,3,1,1)))
forceNetwork(Links = xlink, Nodes = xnode, Source = "source",
             Target = "target", Value = "value", NodeID = "name",
             Nodesize = "size",
             radiusCalculation = "Math.sqrt(d.nodesize)+6",
             Group = "group", opacity = 10, legend = T
             ,zoom=T,opacityNoHove=0)

##################################

g <- graph_from_adjacency_matrix(t(xmat),'undirected')
test <- fastgreedy.community(g)
test <- table(test$membership,xnode[,2])
test2 <- matrix(0,nrow(test),ncol(test),dimnames=dimnames(test))
for(i in 1:nrow(test2)){
  for(j in 1:ncol(test2)){
    m <- colSums(test)[j]
    n <- sum(test) - m
    k <- rowSums(test)[i]
    x <- test[i,j]
    test2[i,j] <- dhyper(x,m,n,k,log=F)
  }
}

##################################

g <- graph_from_adjacency_matrix(t(gmat),'directed')
all_shortest_paths(g,from=V(g)[grep('ERMP1',names(V(g)))],V(g)[grep("signal",names(V(g)))])[[1]]
from <- names(which(gmat[,grep('MIR564',colnames(gmat))]>0))
from <- substr(from,3,nchar(from))
from[from%in%paste(xnode$name)]


x.s <- paste(xnode$group,xnode$name,sep=":")[xlink$source+1]
x.t <- paste(xnode$group,xnode$name,sep=":")[xlink$target+1]
g.s <- paste(xnode$group[xlink$source+1])
g.t <- paste(xnode$group[xlink$target+1])
write.csv(unique(cbind(x.s,x.t,g.s,g.t)),'temp.csv',quote=F,row.names=F)

###############################
#Looking for genes

g <- graph_from_data_frame(gdf)
v_sel <- V(g)[-grep("d'|ph'",names(V(g)))]

len2d <- sapply(v_sel,function(v.f){
  shortest.paths(g,v.f,V(g)[grep("d'",names(V(g)))])
})
path2d <- lapply(v_sel,function(v.f){
  all_shortest_paths(g,v.f,V(g)[grep("d'",names(V(g)))])[[1]]
})

onestep2d <- lapply(1:6,function(x){
  names(which(len2d[x,]==1))
})
names(onestep2d) <- grep("d'",names(V(g)),value=T)
twostep2d <- lapply(onestep2d,function(x){
  v_to <- V(g)[names(V(g))%in%x]
  v_from <- V(g)[-grep("d'|ph'|p'",names(V(g)))]
  rlt <- sapply(v_from,function(v){shortest.paths(g,v,v_to)})
  rlt <- apply(rlt,2,function(x){sum(x==1)})
  test <- rlt[rlt>=quantile(rlt,0.9)]
  # test <- rlt
  test
}) 

library(dplyr)
rlt <- lapply(1:6,function(i){
  x1 <- cbind(onestep2d[[i]],score=99)
  x2 <- cbind(names(twostep2d[[i]]),as.numeric(twostep2d[[i]]))
  out <- as.data.frame(rbind(x1,x2))
  arrange(out %>% group_by(V1) %>% summarise(score=sum(as.numeric(paste(score)))),desc(score))
})
names(rlt) <- names(onestep2d)
rlt <- lapply(rlt,as.data.frame)
write.csv(rlt[[6]],'clipboard',row.names=F)

x <- rlt[[3]]
x1 <- paste(filter(x,score>=99)$V1)
x1 <- paste(x$V1)
x1 <- lapply(c("p'","g'","m'"),function(xi){cbind(gsub("p'|g'|m'","",x1[(grepl(xi,x1))]))})
names(x1) <- c('Pathway','Genotype','Methylation')


######################

g <- graph_from_data_frame(gdf)
from <- V(g)[grep("g'POU",names(V(g)))]
to <- V(g)[grep("d'AD",names(V(g)))]
ps <- all_simple_paths(g,from=from,to=to)
p <- ps[(sapply(ps,function(x){sum(grepl("ph'",names(x)))})<=1)&
          sapply(ps,function(x){sum(grepl("m'",names(x)))})<=1]
p <- p[sapply(p,length)<=5]
# p <- all_shortest_paths(g,from,to)[[1]]
p <- unique(unlist(lapply(p,names)))
p <- t(gmat[rownames(gmat)%in%p,colnames(gmat)%in%p])
plotg(graph_from_adjacency_matrix(p),c(.8,5,1,2))

#######################
#AD Plot

# node.sel <- grepl("d'AD|ph'",colnames(gmat))|colnames(gmat)%in%unique(c(onestep2d[[3]]))
node.sel <- grepl("d'AD|ph'",colnames(gmat))|colnames(gmat)%in%c(onestep2d[[3]],names(which((twostep2d[[3]])>30)))
temp <- gmat[node.sel,node.sel]
g <- graph_from_adjacency_matrix(t(temp))
test <- lapply(paste(rlt[[3]]$V1),function(x){
  print(x)
  from <- V(g)[grep(x,names(V(g)))]
  to <- V(g)[grep("d'AD",names(V(g)))]
  out <- all_simple_paths(g,from=from,to=to)
  out
})
node.sel <- c(names(sort(table(names(unlist(test))),T)[1:10]),
              c(onestep2d[[3]],names(which((twostep2d[[3]])>30))))
node.sel <- colnames(gmat)%in%node.sel

####################################
#Plot pathway network
plotnet <- function(x,mode='undirected',cg=c(.3,1,1,1)){
  diag(x) <- 0
  plot(graph_from_adjacency_matrix(t(x),mode=mode),
       edge.arrow.size=cg[1],
       vertex.size=cg[2],
       vertex.label.cex=cg[3],
       edge.width=cg[4])
}
plotnet(rlt_p2pinp[[2]][[18]],'directed',cg=c(.5,5,.8,1))
