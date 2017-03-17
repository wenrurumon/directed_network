
rm(list=ls())
setwd('C:\\Users\\zhu2\\Documents\\signaling\\codes\\')
source('sparse_2sem_final.R')
source('local_cnif_macro.R')
source('CNIF.R')
sourceCpp("score_function_regression.cpp")
sourceCpp("simple_cycle.cpp")
sourceCpp("initial_sem.cpp")

setwd('C:\\Users\\zhu2\\Documents\\hug')
load('fulldata_adfull.RData')
setwd('C:\\Users\\zhu2\\Documents\\getpathway\\model20170215')
load('C:/Users/zhu2/Documents/getpathway/model20170215/summary/expr_network_20170308.rda')
load("C:/Users/zhu2/Documents/hug/snp_ad.rda")

rlt <- rlt[[115]]

####################
#data process

x <- lapply(snp,function(x){
  x <- x$q2fpca
  x[[1]][,1:which(x[[2]]>=0.9)[1],drop=F]
})
xsn <- sapply(x,ncol)
names(xsn) <- paste('snp',names(xsn),sep='.')
x <- data.matrix(scale(do.call(cbind,x))[,])
rownames(x) <- as.numeric(gsub('[^0-9]','',rownames(snp[[1]][[1]])))

expr_all <- data.matrix(scale(do.call(cbind,lapply(fulldata,function(x) x$expr)))[,])
id.overlap <- rownames(x)[rownames(x)%in%rownames(expr_all)]

x <- x[rownames(x)%in%id.overlap,]
expr_all <- expr_all[rownames(expr_all)%in%id.overlap,]

####################
#train and test

set.seed(1234);test <- sample(1:447,47)
expr_train <- expr_all[-test,]
snp_train <- x[-test,]
expr_test <- expr_all[test,]
snp_test <- x[test,]

set.seed(1234); gc(); system.time(model.sem <- sparse_2sem(Y=expr_train,lambda=0.25,times=100))
gc(); system.time(model.dag <- CNIF(data=expr_train,init.adj=model.sem[[1]]>=0.8,max_parent = 3))

####################
#

gc(); system.time(model.rlt <- sparse_2sem(
  Y=expr_train,Y.fixed=model.dag,
  X=snp_train,lambda=0.3,xsn=xsn,times=100))

plot(graph_from_adjacency_matrix(t(model.dag)))

#####################
#Linear Model Prediction
data_train <- cbind(expr_train,snp_train)
data_test <- cbind(expr_test,snp_test)
pred_lm_test <- sapply(1:ncol(expr_train),function(i){
  print(i)
  yi_train <- data_train[,i]
  xi_train <- data_train[,-i]
  if(ncol(xi_train)==0){
    return(rep(mean(yi_train),nrow(data_test)))
  } else {
    xi_test <- data_test[,-i]
    xi_lm <- coef(lm(yi_train~xi_train))
    as.numeric(data.matrix(cbind(1,xi_test)) %*% matrix(xi_lm,ncol=1))
  }
})

#####################
#DAG Prediction

data_train <- cbind(expr_train,snp_train)
data_test <- cbind(expr_test,snp_test)
pred_dag_test <- sapply(1:ncol(expr_train),function(i){
  print(i)
  yi_train <- data_train[,i]
  xi_train <- data_train[,which(model.rlt$eq_stability[i,]>=0.8),drop=F]
  if(ncol(xi_train)==0){
    return(rep(mean(yi_train),nrow(data_test)))
  } else {
    xi_test <- as.data.frame(cbind(1,data_test[,which(model.rlt$eq_stability[i,]>=0.8),drop=F]))
    xi_lm <- coef(lm(yi_train~xi_train))
    return(as.numeric(data.matrix(xi_test) %*% matrix(xi_lm,ncol=1)))
  }
})

lm_rlt <- cor(pred_lm_test,expr_test)
dag_rlt <- cor(pred_dag_test,expr_test)
rlt2 <- cbind(diag(lm_rlt),diag(dag_rlt))
head(rlt2)

lm_res <- ((pred_lm_test-expr_test)^2) 
dag_res <- ((pred_dag_test - expr_test)^2) 
rlt_res <- (cbind(lm_res=colMeans(lm_res),dag_res=colMeans(dag_res)))
plot(rlt_res[,1])
lines(rlt_res[,2],type='p',col=2)

par(mfrow=c(5,4))
for(j in 1:10){
  plot.ts(rowMeans(expr_test[,j,drop=F]),main='linear model')
  lines(rowMeans(pred_lm_test[,j,drop=F]),col=3)
  plot.ts(rowMeans(expr_test[,j,drop=F]),main='causal network')
  lines(rowMeans(pred_dag_test[,j,drop=F]),col=2)
}

