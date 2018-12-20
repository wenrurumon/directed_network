
rm(list=ls())

setwd('e:\\uthealth\\signaling\\codes\\')
source('sparse_2sem_final.R')
source('local_cnif_macro.R')
source('CNIF.R')
sourceCpp("score_function_regression.cpp")
sourceCpp("simple_cycle.cpp")
sourceCpp("initial_sem.cpp")

setwd('e:\\uthealth\\hug')
load('fulldata_adfull.RData')
setwd('e:\\uthealth\\getpathway\\model20170215')
load('e:\\uthealth\\getpathway/model20170215/summary/expr_network_20170308.rda')
load("e:\\uthealth\\hug/snp_ad.rda")
load("E:/uthealth/getpathway/gene39761.rda")

disease2 <- disease[match(rownames(raw_exp),disease[,1]),-1]
rlt <- rlt[[115]]

##########################
# Macro
##########################

p_ginv_sq <- function(X,p){
  X.eigen = eigen(X);
  X.rank = sum(X.eigen$values>1e-8);
  X.value = X.eigen$values[1:X.rank]^(-1*p);
  if (length(X.value)==1){
    D = as.matrix(X.value);
  }else{
    D = diag(X.value);
  }
  rlt = X.eigen$vectors[,1:X.rank] %*% D %*% t(X.eigen$vectors[,1:X.rank]);
  return(rlt);
}
mrank <- function(X){
  X.svd = svd(X);
  X.rank = sum(X.svd$d>1e-6);
  return(X.rank);
}
mrank_sq <- function(X){
  X.eigen = eigen(X);
  X.rank = sum(Re(X.eigen$values)>1e-6);
  return(X.rank);
}
CCA_chisq_test <- function(rho,n,p,q){
  tstat = -1*n*sum(log(1-rho^2));
  p_value = pchisq(tstat,(p*q),lower.tail=FALSE);
  return(p_value);          
}
cca <- function(A,B){
  n = nrow(A);
  p = mrank(A);
  q = mrank(B);
  if (p <= q){
    X = A;
    Y = B;
  }else{
    X = B;
    Y = A;
  }
  R = p_ginv_sq(cov(X),0.5) %*% cov(X,Y) %*% p_ginv_sq(cov(Y),1) %*% cov(Y,X) %*% p_ginv_sq(cov(X),0.5);
  k = mrank_sq(R);
  d = Re(eigen(R)$values);
  rho = d[1:k]^(0.5);
  rho[rho >= 0.9999]=0.9;
  chisq_p = CCA_chisq_test(rho,n,p,q);
  return(c("rho"=rho,"chisq_p"=chisq_p,"df"=p*q));
}
ccap <- function(A,B){as.numeric(cca(A,B)[2])}
ccaps <- function(As,Bs){
  rlt <- lapply(As,function(A){
    sapply(Bs,function(B){
      ccap(A,B)
    })
  })
  return(unlist(rlt))
}

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
disease2 <- disease[match(rownames(raw_exp),disease[,1]),-1]
disease2[,4] <- as.numeric(disease2[,4]%in%4:5)
d <- disease2[,4]

####################
#train and test

set.seed(1234);test <- sample(1:447,47)
expr_train <- expr_all[-test,]
snp_train <- x[-test,]
expr_test <- expr_all[test,]
snp_test <- x[test,]
d_train <- d[-test]
d_test <- d[test]

set.seed(1234); gc(); system.time(model.sem <- sparse_2sem(Y=expr_train,lambda=0.25,times=100))
gc(); system.time(model.dag <- CNIF(data=expr_train,init.adj=model.sem[[1]]>=0.8,max_parent = 2))

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

#####################
#Validation

lm_rlt <- cor(pred_lm_test,expr_test)
dag_rlt <- cor(pred_dag_test,expr_test)
rlt2 <- cbind(diag(lm_rlt),diag(dag_rlt))
head(rlt2)

lm_res <- (abs(pred_lm_test-expr_test)) 
dag_res <- (abs(pred_dag_test - expr_test)) 
rlt_res <- (cbind(lm_res=colMeans(lm_res),dag_res=colMeans(dag_res)))
plot(rlt_res[,1])
lines(rlt_res[,2],type='p',col=2)

rbind(summary(colMeans(dag_res)),summary(colMeans(lm_res)))
sum(dag_res)/sum(lm_res)

par(mfrow=c(5,4))
for(j in 1:10){
  plot.ts(rowMeans(expr_test[,j,drop=F]),main='linear model')
  lines(rowMeans(pred_lm_test[,j,drop=F]),col=3)
  plot.ts(rowMeans(expr_test[,j,drop=F]),main='causal network')
  lines(rowMeans(pred_dag_test[,j,drop=F]),col=2)
}

#####################################################################
#####################################################################

snp_all <- sapply(snp,function(x){
  x <- x$q2fpca
  x <- x[[1]][match(rownames(expr_all),paste(as.numeric(gsub('[^0-9]','',rownames(snp[[1]][[1]]))))),1,drop=F]
  x
})
colnames(snp_all) <- paste0('snp_',colnames(expr_all))
colnames(expr_all) <- paste0('expr_',colnames(expr_all))
x <- as.data.frame(cbind(d=d,expr_all,snp_all))

library(MASS)
z <- lda(d~.,data=x[-test,])
table(predict(z,x[test,])$class,x$d[test])

######################

test1 <- apply(x[,-1],2,function(i){
  fits(d,i)
})
x2 <- x[,-1][,test1[2,]==0,drop=F]
test2 <- apply(x2[,],2,function(i){fits(d,i)})
i<-2

print(i<<-i+1)
x2 <- x2[,test2[2,]==0,drop=F]
test2 <- apply(x2[-test,],2,function(i){fits(d[-test],i)})
length(test2[2,]);sum(test2[2,]!=0)
z <- lda(d~.,data=as.data.frame(cbind(d,x2[,test2[2,]==0]))[-test,])
table(predict(z,as.data.frame(cbind(d,x2[,test2[2,]==0]))[test,])$class,x$d[test])
