
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

#expression residual
setwd('/home/zhu/rushdata/tacc')
load('expr_residual.rda')
names(expr_residual) <- sapply(expr_residual,colnames)

#df
setwd('/home/zhu/rushdata/tacc/rlt')
load('qtlresult.rda')
df <- cbind(to=paste(rlt[[2]][,1]),from=paste0('m',rlt[[2]][,2]))

#methylation data
setwd('/home/zhu/rushdata/tacc/cut')
files <- dir(pattern='m')
methydata <- lapply(files,function(x){
	load(x)
	return(x)
})
methydata <- do.call(c,methydata)
names(methydata) <- paste0('m',names(methydata))

#Go
edata <- expr_residual[names(expr_residual)%in%df[,1]]
mdata <- methydata[names(methydata)%in%df[,2]]
df <- df[df[,1]%in%names(edata)&df[,2]%in%names(mdata),]
test <- function(to,from){
	d.to <- edata[[which(names(edata)==to)]]
	d.from <- mdata[[which(names(mdata)==from)]]
	ccap2(d.to,d.from)
}
rlt <- sapply(1:nrow(df),function(i){
	print(i)
	x <- df[i,]
	test(x[1],x[2])
})
rlt[is.na(rlt)] <- 1
#rlt2 <- rlt * length(expr_residual) * length(snpdata)
rlt2 <- rlt * nrow(df)

df_m2eres <- df[rlt2<0.01,]
setwd('/home/zhu/rushdata/residual_model')
save(df_m2eres,file='df_m2eres.rda')
