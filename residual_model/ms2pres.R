

rm(list=ls())

#cca

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
ccap <- function(A,B){as.numeric(cca(A,B)[1])}
ccap2 <- function(A,B){
	out <- try(ccap(A,B))
	if(is.numeric(out)){
		return(out)
	} else {
		return(NA)
	}}

#candidate
setwd('/home/zhu/rushdata/residual_model')
for(i in dir()){load(i)}

#snpdata
setwd('/home/zhu/rushdata/tacc/cut')
files <- dir(pattern='snp')
snpdata <- lapply(files,function(x){
	load(x)
	return(x)
})
snpdata <- do.call(c,snpdata)
names(snpdata) <- paste0('s',names(snpdata))

#methylation data
setwd('/home/zhu/rushdata/tacc/cut')
files <- dir(pattern='m')
methydata <- lapply(files,function(x){
	load(x)
	return(x)
})
methydata <- do.call(c,methydata)
names(methydata) <- paste0('m',names(methydata))

#candidate
candidates <- unique(c(df_m2eres[,2],df_s2eres[,2],df_s2mres[,2]))
alldata <- c(snpdata,methydata)
alldata <- alldata[names(alldata)%in%candidates]

#pheres
setwd('/home/zhu/rushdata/tacc')
load('pdata_res.rda')
pid <- rownames(methydata[[1]])
pdata.res <- lapply(pdata.res,function(x){
	x[match(pid,rownames(x)),,drop=F]
})
pheres <- do.call(cbind,pdata.res)

#candidates step1
i <- 0
can1 <- sapply(alldata,function(x){
	print(i<<-i+1)
	ccap2(x,pheres)
})
alldata <- alldata[can1<=0.05]

#step2
i <- 0
test <- sapply(pdata.res,function(x){
	print(i<<-i+1)
	sapply(alldata,function(y){
		ccap2(x,y)
		})
	})
test <- reshape::melt(test)
colnames(test)[1:2] <- c('from','to')
df_ms2pres <- test[test$value < 0.05/length(test$value),]
df_ms2pres <- cbind(to=paste(df_ms2pres[,2]),from=paste(df_ms2pres[,1]))
save(df_ms2pres,file='df_ms2pres.rda')
