rm(list=ls())

##########################
# KEGG mining
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
  unlist(lapply(strsplit(html,": "),function(x){
    out <- try(sapply(strsplit(x[[2]],'\"'),function(x){x[[1]]}))
    if(length(x)>1){return(tolower(out))}
  }))
})
names(pathwaymap) <- tolower(substr(html2,regexpr('">',html2)+2,regexpr('</a>',html2)-1))

##########################
# Pathway Network
##########################

setwd('C:\\Users\\zhu2\\Documents\\getpathway')
load("~/getpathway/gene39761.rda.RData")

enames <- tolower(names(expinpath))[names(expinpath)%in%pathlist[,2]]
expinpath <- expinpath[names(expinpath)%in%pathlist[,2]]

pathwaymap2 <- lapply(unique(c(names(pathwaymap),unlist(pathwaymap))),function(x){
  as.character(unlist(pathwaymap[names(pathwaymap)==x]))
})
names(pathwaymap2) <- unique(c(names(pathwaymap),unlist(pathwaymap)))

test <- sapply(pathwaymap2,function(x){
  unique(c(names(pathwaymap2),unlist(pathwaymap2)))%in%x
})
rownames(test) <- unique(c(names(pathwaymap2),unlist(pathwaymap2)))
tnames <- rownames(test)

enames[!enames%in%tnames]
enames[enames=="glycosylphosphatidylinositol(gpi)-anchor biosynthesis"] <- "glycosylphosphatidylinositol (gpi)-anchor biosynthesis"
enames[enames=="glycosphingolipid biosynthesis - globo series"] <- "glycosphingolipid biosynthesis - globo and isoglobo series"

test <- test[match(enames,rownames(test)),match(enames,colnames(test))]
dimnames(test) <- list(names(expinpath),names(expinpath))
pthnet <- test[!apply(test,1,function(x){all(is.na(x))}),!apply(test,1,function(x){all(is.na(x))})]

##########################
# pathway group network
##########################

pathwaymap <- lapply(pathwaymap2,function(x){unique(x[x%in%enames])})
pathwaymap <- pathwaymap[names(pathwaymap)%in%enames]
pathwaymap <- lapply(pathwaymap,function(x){
  names(expinpath)[match(x,enames)]
})
names(pathwaymap) <- names(expinpath)[match(names(pathwaymap),enames)]

pathlist2 <- pathlist[pathlist[,2]%in%names(expinpath),]

grpmap <- lapply(pathwaymap,function(x){unique(pathlist2[pathlist2[,2]%in%x,1])})
grpname <- pathlist2[match(names(grpmap),pathlist2[,2]),1]
grpmap <- lapply(unique(grpname),function(x){
  unique(unlist(grpmap[grpname==x]))
})
names(grpmap) <- unique(grpname)

test <- do.call(cbind,lapply(grpmap,function(x){
  unique(grpname) %in% x
}))
dimnames(test) <- list(unique(grpname),unique(grpname))
grpnet <- test

save(pathwaymap,pthnet,grpmap,grpnet,file='pathwaymap.rda')
