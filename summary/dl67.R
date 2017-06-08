
###################################################
############## GET DATA ###########################
###################################################

rm(list=ls())
library(data.table)
library(dplyr)
library(igraph)

setwd('C:\\Users\\zhu2\\Documents\\getpathway\\model20170215\\deadline67')
load('genenet.rda'); load('pathnet.rda')
g2d <- read.csv('C:\\Users\\zhu2\\Documents\\getpathway\\model20170215\\gene2disease.csv',stringsAsFactors=F)
g2d <- as.matrix(g2d)
genenet <- as.matrix(genenet)
genenet <- unique(rbind(genenet,g2d))
pathnet <- unique(rbind(pathnet,g2d))

pathnet <- gsub("ph':","ph'",pathnet)
pathnet <- gsub("ph'AGE","ph'Age",pathnet)

node.p <- names(geneincluster)
node.ph <- unique(grep("ph'",as.vector(pathnet),value=T)); node.ph <- substr(node.ph,4,nchar(node.ph))
node.d <- unique(grep("d'",as.vector(pathnet),value=T)); node.d <- substr(node.d,3,nchar(node.d))
pathnet <- gsub("g'","s:",pathnet)
pathnet <- gsub("m'","m:",pathnet)
pathnet <- gsub("d'|ph'|p'","",pathnet)
genenet <- gsub("m::","m:",genenet)
genenet <- gsub("s::","s:",genenet)

gnet <- as.data.table(genenet); colnames(gnet) <- c('from','to')
pnet <- as.data.table(pathnet)
gg <- graph_from_data_frame(gnet)
pg <- graph_from_data_frame(pnet)

###################################################
############## Count GENE #########################
###################################################

#calculate all nodes
node <- (unique(unlist(gnet)))
node.g <- node[!node%in%c(node.d,node.ph)]
node.s <- gsub('s:','',grep('s:',node.g,value=T))
node.m <- gsub('m:','',grep('m:',node.g,value=T))
node.e <- node.g[!grepl('s:|m:',node.g)]
node.g <- table(c(node.s,node.m,node.e))

#node to pheno
gene2d <- filter(gnet,to%in%node.ph & (!from%in%c(node.d,node.ph)))$from
gene2d.s <- gsub('s:','',(unique(grep("s:",gene2d,value=T))))
gene2d.m <- gsub('m:','',(unique(grep("m:",gene2d,value=T))))
gene2d.e <- unique(gene2d[!grepl("m:|s:",gene2d)])
length(gene2d <- table(c(gene2d.s,gene2d.m,gene2d.e)))

#node to expression
gene2d <- filter(gnet,to%in%node.e & grepl('s:|m:',from))
length(unique(gene2d$to))
length(unique(grep('s:',gene2d$from,value=T)))
length(unique(grep('m:',gene2d$from,value=T)))

#node to disease
gene2d <- filter(gnet,to%in%node.d & (!from%in%c(node.d,node.ph)))$from
gene2d.s <- gsub('s:','',(unique(grep("s:",gene2d,value=T))))
gene2d.m <- gsub('m:','',(unique(grep("m:",gene2d,value=T))))
gene2d.e <- unique(gene2d[!grepl("m:|s:",gene2d)])
length(gene2d <- table(c(gene2d.s,gene2d.m,gene2d.e)))

#paths count
gene2d <- filter(gnet,to%in%node.d & (!from%in%c(node.d,node.ph)))

node0 <- unique(grep('s:|m:',filter(gene2d,to%in%node.d)$from,value=T))
cango <- function(g,node1,node2){
  node2 <- node2[!node2%in%node1]
  p <- all_shortest_paths(g,node1,node2)$res
  sapply(p,function(pi) names(pi)[length(pi)])
}
node0.score <- nchar(gsub('[0-9]|m:|s:','',node0))/nchar(gsub('m:|s:','',node0))
node0 <- node0[node0.score>=.66]

#Select node method

g <- pnet
g[which(g$to%in%node.d)]$to <- 'Disease'
node2 <- 'Disease'
imax <- 1
for(i in 1:imax){node2 <- unique(c(node2,filter(g,to%in%node2)$from))}
node1 <- node0
for(i in 1:imax){node1 <- unique(c(node1,filter(g,from%in%node1)$to))}
length(node1); length(node2)
g <- unique(filter(g,from%in%c(node1,node2) & to%in%c(node1,node2)))
g2 <- graph_from_data_frame(g)
node0 <- node0[node0%in%names(V(g2))]

noderm <- c()
for(i in node0){
  print(i)
  if(i %in% noderm){next}
  temp <- unique(cango(g2,i,node0))
  noderm <- c(noderm,temp)
}
node00 <- node0[!node0%in%unlist(noderm)]
nodea <- unique(c(node1,node2))

#Count nodei
cango2 <- function(g,node1,node2){
  length(all_shortest_paths(g,node1,node2)[[1]])
}

w <- sapply(nodea,function(n){
  print(n)
  sapply(nodea,function(n2){
    cango2(g2,n,n2)
  })
})
diag(w) <- 0
sum(w[,colnames(w)%in%node0])

###################################################
############## GET PATH ###########################
###################################################

#Cscore to AD

cscore <- node.ph[c(7,8,12,13,14,15,17,20)]
path_cscore2ad <- lapply(cscore,function(i){all_simple_paths(gg,i,"AD")})

#CREBBP to Cscore

path_crebbp2cscore <- all_shortest_paths(pg,'s:CREBBP',cscore)
# path_crebbp2cscore <- path_crebbp2cscore[sapply(path_crebbp2cscore,length)<=7]
path_crebbp2cscore <- all_shortest_paths(pg,'s:CREBBP',cscore)[[1]]
node.count <- table(unlist(lapply(path_crebbp2cscore,names)))
node.count <- sapply(path_crebbp2cscore,function(x){
  out <- mean(node.count[names(node.count)%in%names(x)])
  names(out) <- names(x)[length(x)]
  out
})
path_crebbp2cscore[paste(names(node.count),node.count)%in%paste(names(tapply(node.count,names(node.count),max)),tapply(node.count,names(node.count),max))]

#POU3F2 to Cscore

# path_crebbp2cscore <- all_shortest_paths(pg,'s:CREBBP',cscore)
# path_crebbp2cscore <- path_crebbp2cscore[sapply(path_crebbp2cscore,length)<=7]
path_pou3f22cscore <- all_shortest_paths(pg,'m:POU3F2',cscore)[[1]]
node.count <- table(unlist(lapply(path_pou3f22cscore,names)))
node.count <- sapply(path_pou3f22cscore,function(x){
  out <- mean(node.count[names(node.count)%in%names(x)])
  names(out) <- names(x)[length(x)]
  out
})
path_pou3f22cscore[paste(names(node.count),node.count)%in%paste(names(tapply(node.count,names(node.count),max)),tapply(node.count,names(node.count),max))]
all_simple_paths(pg,'Olfactory transduction','MMSE')

#Gnet CREBBP
crebbp2 <- strsplit("CRKL CYP3A5 GABRB1 HLA-DRB1 IGSF5 KIR2DL1 MT-ND4 OXSM PIK3CG SQLE",' ')[[1]]
all_shortest_paths(gg,crebbp2[[1]],'Episodic_Memory')[[1]]
