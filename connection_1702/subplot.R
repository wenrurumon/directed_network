
rm(list=ls())
library(igraph)

########################
# SUbnetwork on pathways
########################

fc <- function(x){
  fc <- membership(fastgreedy.community(x))
  fc[] <- match(fc,unique(fc))
  fc
}

setwd('/Users/wenrurumon/Documents/uthealth/subnetplot')
pnet <- read.csv('path_clustering.csv')
load('pmap.rda')
psmap <- pmap
pmap <- unlist(pmap)
pmap <- cbind(group=substr(names(pmap),1,5),path=paste(pmap))
psmap <- cbind(substr(names(psmap),1,5),names(psmap))

g <- graph_from_data_frame(pnet[,1:2],directed=F)
g.fc <- fc(g)

########################
# Plotit
########################

pnode <- unique(
  rbind(cbind(name=paste(pnet[,1]),group=paste(pnet[,3])),
        cbind(paste(pnet[,2]),paste(pnet[,4]))))
p <- gsub(" |,|'",'',pnode[,1])
p2 <- gsub(" |,|'",'',pmap[,2])
group <- psmap[match(pmap[match(p,p2),1],psmap[,1]),2]
group <- substr(group,regexpr(' ',group)+1,nchar(group))
pnode[,2] <- group

pnode <- data.frame(pnode,size=5)
plink <- data.frame(source=match(pnet[,1],pnode$name)-1,
                    target=match(pnet[,2],pnode$name)-1,value=1)
library(networkD3)
theplot <- forceNetwork(Links = plink, Nodes = pnode, Source = "source",
             Target = "target", Value = "value", NodeID = "name",
             Nodesize = "size",linkColour = "#999",
             radiusCalculation = "Math.sqrt(d.nodesize)+6",
             Group = "group", opacity =1,charge=-5, legend = T
             ,zoom=T,opacityNoHove=1)
theplot
#saveNetwork(theplot,'pathwayclustering.html')
