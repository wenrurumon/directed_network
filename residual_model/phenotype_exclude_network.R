
rm(list=ls())
load("C:/Users/zhu2/Documents/getpathway/model20170215/summary/phenet.rda")

data <- phenet$data
adj <- phenet$adj
colnames(adj) <- rownames(adj)
adj <- dplyr::filter(reshape::melt(adj),value==1)
colnames(adj)[1:2] <- c('to','from')

themodel <- function(iy,ix,data){
  data <- data.matrix(data)
  y <- data[,iy]
  if(length(ix)==0){
    res <- y
  } else {
    x <- data[,ix,drop=F]
    res <- lm(y~x)$residual
  }
  res <- cbind(res)
  colnames(res) <- colnames(data)[iy]
  res
}

p <- colnames(data)
pdata.res <- lapply(p,function(pi){
  iy <- which(p==pi)
  px <- adj[which(adj[,1]==pi),2]
  ix <- which(p%in%px)
  themodel(iy,ix,data)
})
names(pdata.res) <- sapply(pdata.res,colnames)
save(pdata.res,file='pdata_res.rda')
