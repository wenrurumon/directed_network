
rm(list=ls())
library(igraph)
library(dplyr)
library(reshape)

mat2df <- function(mat){
  rlt <- filter(melt(mat),value>0)[,1:2,drop=F]
  colnames(rlt) <- c('to','from')
  rlt
}
path2df <- function(pathi){
  v <- names(pathi)
  rlt <- do.call(rbind,lapply(1:(length(v)-1),function(i){v[i:(i+1)]}))
  colnames(rlt) <- c('from','to')
  rlt
}

setwd('C:\\Users\\zhu2\\Documents\\getpathway\\model20170215\\tacc')
load(file='expr_residual.rda')
load('C:/Users/zhu2/Documents/getpathway/model20170215/rlt_cluster_network.rda')
rlt.g2g <- list(data=exprincluster,adj=rlt_all_split,adj_c2c=cluster_network)
rlt.g2g$df <- (do.call(rbind,lapply(rlt.g2g$adj,mat2df)))
rlt.g2g$df <- cbind(to = paste(rlt.g2g$df[,1]), from = paste( rlt.g2g$df[,2]))
rlt.g2g$df <- rlt.g2g$df[rlt.g2g$df[,1]!=rlt.g2g$df[,2],,drop=F]
load("C:/Users/zhu2/Documents/getpathway/model20170215/summary_201705/rltsel.rda")
geneincluster <- lapply(exprincluster,colnames)

#gene to pathway

g <- igraph::graph_from_data_frame(data.frame(rltsel[,2],rltsel[,1]))
to <- unique(unlist(geneincluster[grep('AMPK',names(geneincluster))]));to
to <- to[to%in%names(V(g))]
from <- 'm:KIF4B'
from <- grep(from,names(V(g)),value=T); from
paths <- all_shortest_paths(g,from[1],to)[[1]]
paths <- paths[sapply(paths,length)<=median(sapply(paths,length))]
gene.next <- unique(names(sapply(paths,function(x){x[length(x)]})))

#pathway to another pathway

from <- paths[[1]]
to <- unique(unlist(geneincluster[grep('PI3K',names(geneincluster))]))
to <- to[to%in%names(V(g))]
paths <- all_simple_paths(g,from[1],to)[[1]]
paths <- paths[sapply(paths,length)==max(sapply(paths,length))]
paths <- unique(names(sapply(paths,function(x){x[length(x)]})))

#To disease

from <- paths[[1]]
to <- 'AD'
paths <- all_simple_paths(g,from[1],to)[[1]]

##################################################################################
##################################################################################
##################################################################################
#POU3F2 to AD
rm(list=ls())
library(igraph)
library(dplyr)
library(reshape)
mat2df <- function(mat){
  rlt <- filter(melt(mat),value>0)[,1:2,drop=F]
  colnames(rlt) <- c('to','from')
  rlt
}
path2df <- function(pathi){
  v <- names(pathi)
  rlt <- do.call(rbind,lapply(1:(length(v)-1),function(i){v[i:(i+1)]}))
  colnames(rlt) <- c('from','to')
  rlt
}
checkpathway <- function(x,y=NULL){
  x <- names(geneincluster)[sapply(geneincluster,function(y){x%in%y})]
  if(is.null(y)){return(x)}
  grep(y,x,value=T)
}
setwd('C:\\Users\\zhu2\\Documents\\getpathway\\model20170215\\tacc')
load(file='expr_residual.rda')
load('C:/Users/zhu2/Documents/getpathway/model20170215/rlt_cluster_network.rda')
rlt.g2g <- list(data=exprincluster,adj=rlt_all_split,adj_c2c=cluster_network)
rlt.g2g$df <- (do.call(rbind,lapply(rlt.g2g$adj,mat2df)))
rlt.g2g$df <- cbind(to = paste(rlt.g2g$df[,1]), from = paste( rlt.g2g$df[,2]))
rlt.g2g$df <- rlt.g2g$df[rlt.g2g$df[,1]!=rlt.g2g$df[,2],,drop=F]
load("C:/Users/zhu2/Documents/getpathway/model20170215/summary_201705/rltsel.rda")
geneincluster <- lapply(exprincluster,colnames)
gdf <- data.frame(from=rltsel[,2],to=rltsel[,1],stringsAsFactors=F)
# filter(gdf,to=='AD')


g <- igraph::graph_from_data_frame(data.frame(rltsel[,2],rltsel[,1]))
#POU3F2 to AMPK
to <- unique(unlist(geneincluster[grep('AMPK',names(geneincluster))]))
paths <- all_shortest_paths(g,from='m:KIF4B',to)[[1]]
paths <- paths[sapply(paths,length)<8]
nnext1 <- nnext <- unique(names(sapply(paths,function(x){x[length(x)]})))
passby <- sort(table(names(unlist(paths))))
passby <- sapply(paths,function(pathi){
  min(passby[match(names(pathi),names(passby))])  
})
# tapply(passby,names(sapply(paths,function(x){x[length(x)]})),min)
passby <- unique(unlist(lapply(paths,names)))
#node out 3
genes <- c(passby,unique(unlist(geneincluster[grep('AMPK|PI3K',names(geneincluster))])))
genes <- c(genes,'HDLratio','CHL','LDL','HDL','AD')
g2 <- dplyr::filter(gdf,from%in%genes & to%in%genes)
g2 <- graph_from_data_frame(g2)
all_shortest_paths(g2,'PIK3CA',c('LDL','HDL','CHL','HDLratio','AD'))
all_shortest_paths(g,'s:POU3F2','CREB1')[[1]]

#############################
#############################
#CREBBP
#############################

test <- strsplit("CRKL CYP3A5 GABRB1 HLA-DRB1 IGSF5 KIR2DL1 MT-ND4 OXSM PIK3CG SQLE",' ')[[1]]
test <- lapply(test,function(x){
  rlt <- list()
  rlt[[1]] <- checkpathway(x)
  # rlt[[1]] <- grep('MAPK|FoxO|Jak|PI3K|Hippo',checkpathway(x),value=T)
  rlt[[2]] <- all_shortest_paths(g,x,'AD')[[1]]
  rlt[[3]] <- '###################'
  rlt
})
names(test) <- strsplit("CRKL CYP3A5 GABRB1 HLA-DRB1 IGSF5 KIR2DL1 MT-ND4 OXSM PIK3CG SQLE",' ')[[1]]
allcrebbpqtl <- strsplit("ADRA1B APOL1 ARG1 BPGM CACNG7 CALM2 CAMP CCL5 CFB CHMP4B CLDN11 CLIP1 CREB3 CYP46A1 DNAH14 DNMT1 EIF3H FCF1 FFAR2 GABRB1 GCH1 GRIK5 GUCY2D IFI30 IFNAR2 IRF1 KCNMA1 MIR181D MRC2 OR2M4 PDGFA PIK3CG PRMT5 PSMB9 PSME2 PTMA SEPSECS SMAD9 SNRPE SPOPL SQRDL STX2 SUMO2 TAP1 TAPBP TBXA2R TEK THOC7 XRN1"," ")

passby <- c('CRKL','MT-ND4','OXSM','PIK3CG')
gene1 <- unique(unlist(geneincluster[grep('FoxO',names(geneincluster))]))
gene2 <- unique(unlist(geneincluster[grep('Hippo|Jak-STAT',names(geneincluster))]))
pass1 <- lapply(passby,function(f){
  x <- all_shortest_paths(g,f,gene1)[[1]]
  l <- sapply(x,length)
  if(length(l)==0){return(NULL)}else{return(x[l==min(l)])}
})
nnext <- unique(names(unlist(do.call(c,pass1))))
pass2 <- lapply(nnext,function(f){
  x <- all_shortest_paths(g,f,gene2)[[1]]
  l <- sapply(x,length)
  if(length(l)==0){return(NULL)}else{return(x[l==min(l)])}
})
nnext2 <- unique(names(unlist(do.call(c,pass2))))
d <- c('MMSE','Working_Memory','CHL','Weight','AD')
gene3 <- unique(c(nnext,nnext2))
pass3 <- lapply(gene3,function(f){
  print(f)
  x <- all_shortest_paths(g,f,d)
  x[[1]]
})
geneall <- unique(names(unlist(do.call(c,pass3))))
qtl <- allcrebbpqtl[[1]][which(allcrebbpqtl[[1]]%in%geneall)]
gdf <- data.frame(from=rltsel[,2],to=rltsel[,1],stringsAsFactors=F)
geneall <- c(nnext,nnext2,d)
g2 <- filter(gdf,from%in%geneall&to%in%geneall)
g2 <- rbind(g2,cbind(from="s'CREBBP",to=passby))
g2 <- graph_from_data_frame(g2)
#CREBBP Go

#from 0 to 1
gene021 <- lapply(passby,function(f){
  x <- all_shortest_paths(g2,f,gene1[gene1%in%names(V(g2))])
  x[[1]]
})
pass021 <- do.call(c,gene021)
gene021 <- unique(names(sapply(pass021,function(x){x[length(x)]})))
gene021_path <- lapply(gene021,checkpathway,y='FoxO|Hippo|Jak|Alzheimer')
gene021 <- gene021[sapply(gene021_path,length)>1]
pass021 <- lapply(passby,function(f){
  all_shortest_paths(g2,f,gene021)[[1]]
})

#from 1 to 2
pass122 <- shortest.paths(g2,gene021,gene2[gene2%in%names(V(g2))])
gene122 <- colnames(pass122)[apply(pass122,2,min)<=1]
gene122_path <- lapply(gene122,checkpathway,y='Hippo|Jak|Alzheimer')
gene122 <- gene122[sapply(gene122_path,length)>1]

#from 2 to d
pass22d <- shortest.paths(g2,gene122,d)
gene122_path[sapply(gene122_path,length)>1]
