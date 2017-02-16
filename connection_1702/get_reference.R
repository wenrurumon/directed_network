rm(list=ls())

##########################
# Pathway Group
##########################

setwd("C:/Users/zhu2/Documents/getpathway/model20170215/")

url2 <- 'http://www.genome.jp/kegg/pathway.html'
html2 <- readLines(url2)
# html2 <- grep('map=hsa|map=map',html2,value=T)
html2 <- grep('/kegg-bin/show_pathway?',html2,value=T)

url2 <- paste0('http://www.genome.jp/',
               substr(html2,regexpr("/kegg-bin",html2),regexpr('&amp|">',html2)-1)
)

pathwaymap <- lapply(url2,function(url){
  print(url)
  html <- readLines(url)
  html <- html[which(grepl('show_pathway',html)&grepl('rect',html)&grepl('hsa|map',html))]
  sapply(strsplit(html,": "),function(x){
    out <- try(sapply(strsplit(x[[2]],'\"'),function(x){x[[1]]}))
    ifelse(length(x)==1,x,out)
  })
})
names(pathwaymap) <- substr(html2,regexpr('">',html2)+2,regexpr('</a>',html2)-1)

save(pathwaymap,file="pathwaymap.rda")

##########################
# Pathway Group
##########################

rm(list=ls())

setwd('C:\\Users\\zhu2\\Documents\\getpathway')
load("~/getpathway/gene39761.rda.RData")
load("C:/Users/zhu2/Documents/getpathway/ptwmap2.rda")

setwd("C:/Users/zhu2/Documents/getpathway/model20170215/")
# load('pathwaymap.rda')

pathlist <- pathlist[names(expinpath)%in%pathlist[,2],]
pathwaymap <- lapply(pathwaymap,function(x){
  x[x%in%names(expinpath)]
})[names(pathwaymap)%in%names(expinpath)]

ref.grpnet <- lapply(pathwaymap,function(x){
  unique(pathlist[pathlist[,2]%in%x,1])
})
names(ref.grpnet) <- pathlist[match(names(pathwaymap),pathlist[,2]),1]
refnames <- unique(names(ref.grpnet))[!is.na(unique(names(ref.grpnet)))]

ref.grpnet <- lapply(unique(names(ref.grpnet)),function(x){
  unique(unlist(ref.grpnet[names(ref.grpnet)==x]))
})[!is.na(unique(names(ref.grpnet)))]
names(ref.grpnet) <- refnames

grpnet <- matrix(0,length(ref.grpnet),length(ref.grpnet),dimnames=list(refnames,refnames))
for(i in 1:length(ref.grpnet)){
  grpnet[i,refnames %in% ref.grpnet[[i]]] <- 1
}
# 
# setwd('C:\\Users\\zhu2\\Documents\\getpathway')
# load("~/getpathway/gene39761.rda.RData")
# load("C:/Users/zhu2/Documents/getpathway/ptwmap2.rda")
# 
# expinpath <- expinpath[-3]
# inputs <- lapply(unique(pathlist[,1]),function(grpi){
#   paths <- pathlist[pathlist[,1]==grpi,2]
#   input <- expinpath[names(expinpath)%in%paths]
#   print(list(groupname=grpi,paths=names(input)))
#   input
# })
# names(inputs) <- unique(pathlist[,1])
# inputs <- inputs[sapply(inputs,length)>0]

names1 <- unlist(lapply(inputs,names))
pthnet <- matrix(0,length(names1),length(names1),dimnames=list(names1,names1))
for(i in 1:length(pathwaymap)){
  pthnet[i,names1%in%pathwaymap[[i]]] <- 1
}

###########################
# Summary
###########################

plotnet <- function(x,mode='directed'){
  diag(x) <- 0
  plot(graph_from_adjacency_matrix(t(x),mode=mode),
       edge.arrow.size=.3,
       vertex.size=3,
       vertex.label.cex=1,
       edge.width=.1)
}

plotnet(grpnet)
save(pathwaymap,ref.grpnet,pthnet,grpnet,file="pathwaymap.rda")

