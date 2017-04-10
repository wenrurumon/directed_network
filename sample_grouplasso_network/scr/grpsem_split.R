
rm(list=ls())

setwd('/home/zhu/rushdata/expression_clustering')
library(Rcpp)
library(grplasso)
library(data.table)

load('expr_cluster.rda')

args <- as.numeric(commandArgs(trailingOnly=TRUE))
if(length(args)==0){args<-1}

##############################
# Macro
##############################

#expand variable x to equation set
expandx <- function(x,n){
  temp <- lapply(1:n,function(i){
    temp <- matrix(0,nrow=nrow(x)*n,ncol=ncol(x))
    temp[1:nrow(x)+(i-1)*nrow(x),] <- as.matrix(x)
    temp
  })
  do.call(cbind,temp)
}
#prepare data for equation set
data4sem <- function(y,x){
  x2 <- lapply(x,function(xi){
    expandx(xi,ncol(y))
  })
  list(
    y=as.vector(y),
    x=do.call(cbind,x2),
    index=rep(1:length(x2),sapply(x2,ncol))
  )
}
#quick group lasso sturcture
quick_grplasso <- function(y,x,index,lambda){
  lambda <- lambdamax(x=as.matrix(x),y=as.matrix(y),index=index, 
                      penscale = sqrt, model = LinReg(),
                      center=TRUE,standardized=TRUE) * 0.5^(1/lambda-1)
  fit <- grplasso(x=as.matrix(x),y=as.matrix(y),index=index,lambda=lambda,model=LinReg(),
                  penscale = sqrt,
                  control = grpl.control(update.hess = "lambda", trace = 0))
  tapply(coef(fit),index,function(x){any(x!=0)})
}

###########################################################
#raw group lasso sem within Ys
###########################################################
Y <- lapply(expr_cluster,function(x){x[,1:min(20,ncol(x)),drop=F]})
#system.time(test <- grpsem_split(Y,lambda=0.7))
grpsem_split <- function(Y,lambda=.5,outaslist=FALSE,start=1,end=length(Y)){
  models <- lapply(start:end,function(i){
    print(i)
    data <- data4sem(y=Y[[i]],x=Y[-i])
    coef <- rep(0,length(Y))
    coef[-i] <- quick_grplasso(y=data$y,x=cbind(1,data$x),index=c(NA,data$index),lambda)
    gc();
    as.numeric(coef)
  })
  if(!outaslist){
    return(do.call(rbind,models))
  } else {
    return(models)
  }
}

trail <- seq(1,663,35)
trail[length(trail)] <- 663

test <- grpsem_split(Y,lambda=0.9,FALSE,trail[args],trail[args+1])
save(test,file=paste0('rlt_sem_',args,'.rda'))