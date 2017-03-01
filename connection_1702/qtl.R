rm(list=ls())

library(igraph)
library(WGCNA)
library(grplasso)

##########################
# Macro
##########################

qpca <- function(A,rank=0,ifscale=TRUE){
  if(ifscale){A <- scale(as.matrix(A))[,]}
  A.svd <- svd(A)
  if(rank==0){
    d <- A.svd$d
  } else {
    d <- A.svd$d-A.svd$d[min(rank+1,nrow(A),ncol(A))]
  }
  d <- d[d > 1e-8]
  r <- length(d)
  prop <- d^2; info <- sum(prop)/sum(A.svd$d^2);prop <- cumsum(prop/sum(prop))
  d <- diag(d,length(d),length(d))
  u <- A.svd$u[,1:r,drop=F]
  v <- A.svd$v[,1:r,drop=F]
  x <- u%*%sqrt(d)
  y <- sqrt(d)%*%t(v)
  z <- x %*% y
  rlt <- list(rank=r,X=x,Y=y,Z=x%*%y,prop=prop,info=info)
  return(rlt)
}

qq_plot <- function(p_value){
  n = length(p_value);
  exp = -log10((c(1:n)-0.5)/n);
  rgen = -log10(sort(p_value));
  plot(exp,rgen,xlab="-log10(Expect)",ylab="-log10(Real)");
  abline(0,1,col="red")
}

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
  return(c("chisq_p"=chisq_p,"df"=p*q));
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

######################################
# Load data
######################################

setwd('C:\\Users\\zhu2\\Documents\\getpathway\\model20170215')
load('rlt_p2pinp.rda')
load('disease.rda')
disease[,5] <- ifelse(disease[,5]%in%c(4,5),1,0); disease <- disease[,-1]
load('pid.rda')
setwd('sample_methy_snp')
for(i in dir(pattern='rda')){load(i)}

pca.pathways <- do.call(c,rlt_p2pinp[[1]])
names(pca.pathways) <- do.call(c,lapply(rlt_p2pinp[[1]],names))

y <- do.call(cbind,pca.pathways)
y.qpca <- qpca(y)

####################################
# after loading
####################################

i <- 0
paths <- lapply(pca.pathways,function(x){
  x <- x[,]
  yx <- y.qpca$X[,1:which(y.qpca$prop>.5)[1],drop=F]
  as.matrix(lm(x~yx)$residual)[,,drop=F]
})
cca_d_p0 <- do.call(rbind,lapply(lapply(1:5,function(i){disease[,i,drop=F]}),function(x){
  i <<- i+1
  j <- 0
  t(sapply(paths,function(y){
    j <<- j+1
    c(d=i,p=j,cca(x,y))
  }))
}))

i <- 0
paths <- lapply(pca.pathways,function(x){
  x <- x[,]
  yx <- y.qpca$X[,1:which(y.qpca$prop>.2)[1],drop=F]
  as.matrix(lm(x~yx)$residual)[,,drop=F]
})
cca_d_p1 <- do.call(rbind,lapply(lapply(1:5,function(i){disease[,i,drop=F]}),function(x){
  i <<- i+1
  j <- 0
  t(sapply(paths,function(y){
    j <<- j+1
    c(d=i,p=j,cca(x,y))
  }))
}))

i <- 0
paths <- lapply(pca.pathways,function(x){
  x <- x[,]
  yx <- y.qpca$X[,1:which(y.qpca$prop>.3)[1],drop=F]
  as.matrix(lm(x~yx)$residual)[,,drop=F]
})
cca_d_p2 <- do.call(rbind,lapply(lapply(1:5,function(i){disease[,i,drop=F]}),function(x){
  i <<- i+1
  j <- 0
  t(sapply(paths,function(y){
    j <<- j+1
    c(d=i,p=j,cca(x,y))
  }))
}))

i <- 0
paths <- lapply(pca.pathways,function(x){
  x <- x[,]
  yx <- y.qpca$X[,1:which(y.qpca$prop>.4)[1],drop=F]
  as.matrix(lm(x~yx)$residual)[,,drop=F]
})
cca_d_p3 <- do.call(rbind,lapply(lapply(1:5,function(i){disease[,i,drop=F]}),function(x){
  i <<- i+1
  j <- 0
  t(sapply(paths,function(y){
    j <<- j+1
    c(d=i,p=j,cca(x,y))
  }))
}))

cca_d_p <- list(cca_d_p0,cca_d_p1,cca_d_p2,cca_d_p3)
par(mfrow=c(2,2))
for(i in cca_d_p){
  qq_plot(i[,3])
}
par(mfrow=c(1,1))

cca_d_p0[cca_d_p0[,3] < (0.05/nrow(cca_d_p0)),]

####################################

qtl <- function(i,data=pmethy,maxprop=0.95,pathexclude=0){
  # i <- 1; data <- psnp; maxprop <- .95; pathexclude <- 0.3
  if(pathexclude==0){
    paths <- pca.pathways
  } else {
    paths <- lapply(pca.pathways,function(x){
      x <- x[,]
      yx <- y.qpca$X[,1:which(y.qpca$prop>pathexclude)[1],drop=F]
      as.matrix(lm(x~yx)$residual)[,,drop=F]
    })
  }
  xi <- data[[i]]
  if(is.null(xi)){return(c(p=NA,g=i,rep(NA,3)))}
  pid2 <- rownames(xi[[1]])
  pid2 <- as.numeric(gsub('ROS|MAP','',pid2))
  fpca <- xi$fpca$score[match(pid,pid2),1:which(xi$fpca$prop>maxprop)[1],drop=F]
  qfpca <- xi$qfpca$score[match(pid,pid2),1:which(xi$qfpca$prop>maxprop)[1],drop=F]
  q2fpca <- xi$q2fpca$score[match(pid,pid2),1:which(xi$q2fpca$prop>maxprop)[1],drop=F]
  j <- 0
  rlt <- do.call(rbind,lapply(paths,function(pcai){
    pcai <- pcai[,,drop=F]
    # print(j <<- j+1)
    rbind(c(p=j,g=i,m=1,cca(pcai,fpca)),
          c(p=j,g=i,m=2,cca(pcai,qfpca)),
          c(p=j,g=i,m=3,cca(pcai,q2fpca)))
  }))
  rlt
}

# test <- qtl(42,data=psnp)
datai <- pmethy[1:100]
system.time(test1 <- lapply(1:length(datai),function(i){print(i);qtl(i,data=datai,maxprop=.9,pathexclude=0.3)}))
system.time(test2 <- lapply(1:length(datai),function(i){print(i);qtl(i,data=datai,maxprop=.9,pathexclude=0.2)}))
test <- do.call(rbind,test2)[,4]; qq_plot(test[!is.na(test)])
system.time(test3 <- lapply(1:length(datai),function(i){print(i);qtl(i,data=datai,maxprop=.9,pathexclude=0.4)}))
test <- do.call(rbind,test3); qq_plot(test[!is.na(test)])

#snp with 0.2




