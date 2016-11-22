
#####################################
#Setup
#####################################

library(fda)
library(MASS)
library(flare)
library(corpcor) ###fast.svd only use when data is very large
### Y and X are the data used to construct structure equation model, 
### We use ADMM to solve sparsity of each structure equation, 
### lambda is the parameter for sparse penalty, here we use lasso as penalty function,
### rho is the penalty parameter used in ADMM
### xsn is the parameter to specify the number of scores within each gene region, 
### it can be a vector and its length equals to the  number of genes,
### if it is a vector without names, the default names for each gene is x1 to x-genenumber
### the default value of xsn is 1, which means each variable of X represent a gene,
### stability is the threshold set for stability selection, the default value is 0.8,
### times is the resample times when use stability selection,
### only when set times larger than 1 the preogram will do stability selection.
### the final structure is based on stability.
sparse_2sem<-function(Y,X=NULL,Y.fixed=NULL,XY.fixed=NULL,lambda=0.1,rho=1,xsn=1,stability=0.8,times=1){
  
  #   Y = raw
  #   X <- Y.fixed <- XY.fixed <- NULL
  #   lambda = 0.1
  #   rho = 1
  #   xsn = 1
  #   stability = 0.8
  #   times = 1
  
  Y<-as.matrix(Y)
  n<-nrow(Y)
  m<-ncol(Y)
  if(is.null(X)){mm=0} else {mm<-ncol(X)}
  M<-m+mm
  if(times==1){
    rlt_coef<-c()
    rlt_rs<-c()
    if(is.null(X)){
      ###eq_matrix saves the selected predictors for each equation
      eq_matrix<-matrix(0,m,M)
      if(!is.null(Y.fixed)){
        eq_matrix[,1:m]<-Y.fixed
      }else{
        for(k in 1:m){
          dataX=as.matrix(Y[,-k])
          #xy0=t(dataX)%*%Y[,k]
          #idx=sort.int(xy0,decreasing=TRUE,index.return=TRUE)$ix
          #dataX=dataX[,idx]
          outl1=slim(X=dataX,Y=Y[,k],lambda,rho,method = "lasso",verbose=FALSE)
          ##save data for linear model
          x_sel<-which(outl1$beta!=0)
          x_num<-x_sel+as.numeric(x_sel>=k)
          eq_matrix[k,x_num]<-1
        }
      }
      colnames(eq_matrix)<-colnames(Y)
      rownames(eq_matrix)<-colnames(Y) 
      ##construct equations
      for(k in 1:m){
        rlt_sp<-c()
        x_sel<-which(eq_matrix[k,]!=0)
        if(length(x_sel)>0){
          W=as.matrix(Y[,which(eq_matrix[k,]!=0)])
          coef=ginv(t(W)%*%W)%*%t(W)%*%Y[,k]
          sigma<-(1/n)*sum((Y[,k]-W%*%coef)^2)
          SIGMA<-sigma*ginv(t(W)%*%W)
          for(l in 1:length(x_sel)){
            Tc=coef[l]^2/diag(SIGMA)[l]
            praw<-pchisq(Tc,df=1,lower.tail=FALSE)
            #rlt_sp<-c(rlt_sp,min(praw*length(x_sel),1))##p-value for each score or variable selected in this equation
            rlt_sp<-c(rlt_sp,praw)
          }
          rlt_coef<-rbind(rlt_coef,cbind(coef,rlt_sp))
          rlt_rs<-c(rlt_rs,sigma)
          #eq_system<-paste("Y[,",k,"]~as.matrix(Y[,which(eq_matrix[",k,",]!=0)])",sep="")
          #lm.fit<-lm(eval(parse(text=eq_system)))
          #rlt_coef<-rbind(rlt_coef,summary(lm.fit)$coefficients[-1,])
          #rlt_rs<-c(rlt_rs,(1/n)*sum(summary(lm.fit)$residuals^2))
        }
        else {k<-k+1}
      }
    }
    else{
      rlt_p<-c()##save p-values for each gene calculated from t2 test
      rlt_sp<-c()##save p-values for each score calculated from t2 test
      ###Structural equation model 
      ##first stage
      Yhat<-matrix(0,n,m)
      #Xinv<-ginv(t(X)%*%X)
      Xinv<-myinv(t(X)%*%X)
      Yhat<-X%*%Xinv%*%t(X)%*%Y
      colnames(Yhat)<-colnames(Y)
      ##second stage
      ###eq_matrix saves the selected predictors for each equation
      eq_matrix<-matrix(0,m,M)
      ##beta_matrix saves the coefficients from lasso on each equation
      beta_matrix<-matrix(0,m,M-1)
      for(k in 1:m){
        dataX=cbind(Yhat[,-k],X)
        dataX<-as.matrix(dataX)
        varx<-apply(dataX,2,var)
        outl1=slim(X=dataX[,varx!=0],Y=Y[,k],lambda,rho,method = "lasso",verbose=FALSE)
        
        ##save data for linear model or linear equation
        beta_matrix[k,varx!=0]<-outl1$beta
        x_sel<-which(beta_matrix[k,]!=0)
        x_num<-x_sel+as.numeric(x_sel>=k)
        eq_matrix[k,x_num]<-1
      }
      if(!is.null(Y.fixed)){
        eq_matrix[,1:m]<-Y.fixed
      }
      if(!is.null(XY.fixed)){
        eq_matrix[,(m+1):M]<-XY.fixed
      }
      ##construct equations
      YX<-cbind(Y,X)
      colnames(eq_matrix)<-colnames(YX)
      rownames(eq_matrix)<-colnames(Y) 
      
      for(k in 1:m){
        x_sel<-which(eq_matrix[k,]!=0)
        if(length(x_sel)>0){
          W<-YX[,x_sel]
          What<-as.matrix(X%*%Xinv%*%t(X)%*%W)
          WX<-as.matrix(t(W)%*%X%*%Xinv%*%t(X))
          coef<-ginv(WX%*%W)%*%WX%*%Y[,k]
          sigma<-(1/n)*sum((Y[,k]-W%*%coef)^2)
          SIGMA<-sigma*ginv(t(What)%*%What)
          rlt_coef<-c(rlt_coef,coef)
          rlt_rs<-c(rlt_rs,sigma)
          ##test using t2 statistics for each variable selected in this equation
          pTc<-c()
          for(l in 1:length(x_sel)){
            Tc=coef[l]^2/diag(SIGMA)[l]
            praw<-pchisq(Tc,df=1,lower.tail=FALSE)
            #rlt_sp<-c(rlt_sp,min(praw*length(x_sel),1))##p-value for each score or variable selected in this equation
            rlt_sp<-c(rlt_sp,praw)
          }
          if(sum(xsn)==1){
            rlt_p<-rlt_sp
          }else{##for cases that one gene have more than one scores
            if(is.null(names(xsn))){
              names(xsn)<-paste("x",1:length(xsn),sep="")
            }
            xsn1<-c(rep(1,m),xsn)
            xsn_num<-cumsum(xsn1)
            g_sel<-intersect(which(0<x_sel),which(x_sel<=xsn_num[1]))
            if(length(g_sel)==1){
              Tg=coef[g_sel]^2/diag(SIGMA)[g_sel]
              praw<-pchisq(Tg,df=1,lower.tail=FALSE)
              #rlt_p<-c(rlt_p,min(praw*length(x_sel),1))
              rlt_p<-c(rlt_p,praw)
            }
            if(length(g_sel)>1){
              Tg=t(coef[g_sel])%*%ginv(SIGMA[g_sel,g_sel])%*%coef[g_sel]
              praw<-pchisq(Tg,df=length(g_sel),lower.tail=FALSE)
              #rlt_p<-c(rlt_p,min(praw*length(x_sel),1))
              rlt_p<-c(rlt_p,praw)
            }
            for(j in 2:length(xsn1)){
              g_sel<-intersect(which(xsn_num[j-1]<x_sel),which(x_sel<=xsn_num[j]))
              if(length(g_sel)==1){
                Tg=coef[g_sel]^2/diag(SIGMA)[g_sel]
                praw<-pchisq(Tg,df=1,lower.tail=FALSE)
                #rlt_p<-c(rlt_p,min(praw*length(x_sel),1))
                rlt_p<-c(rlt_p,praw)
              }
              if(length(g_sel)>1){
                Tg=t(coef[g_sel])%*%ginv(SIGMA[g_sel,g_sel])%*%coef[g_sel]
                praw<-pchisq(Tg,df=length(g_sel),lower.tail=FALSE)
                #rlt_p<-c(rlt_p,min(praw*length(x_sel),1))
                rlt_p<-c(rlt_p,praw)
              }
            }
          }
        }
        else {k<=k+1}
      }
      
    }
    ####save gene_score network edge
    varnames<-colnames(eq_matrix)
    Snet<-c()
    for(i in 1:nrow(eq_matrix)){
      for(j in 1:ncol(eq_matrix)){
        if(eq_matrix[i,j]!=0){
          Snet<-rbind(Snet,c(i,varnames[i],j,varnames[j]))
        }
      }
    }
    colnames(Snet)<-c("responseNum","response","predictorNum","predictor")
    ###based on eq_matrix, construct gene network edge
    Gnet<-c()
    if(sum(xsn)==1){
      Gnet<-Snet
    }
    else{
      nx<-length(xsn)
      G_matrix<-matrix(0,m,m+nx)
      G_matrix[,1:m]<-eq_matrix[,1:m]
      score_num<-m
      for(i in 1:nx){
        G_score<-as.matrix(eq_matrix[,(score_num+1):(score_num+xsn[i])])
        G_matrix[,m+i]<-apply(G_score,1,max)
        score_num<-score_num+xsn[i]
      }
      colnames(G_matrix)<-c(colnames(Y),names(xsn))
      rownames(G_matrix)<-colnames(Y)
      Gvarnames<-colnames(G_matrix)
      for(i in 1:nrow(G_matrix)){
        for(j in 1:ncol(G_matrix)){
          if(G_matrix[i,j]!=0){
            Gnet<-rbind(Gnet,c(i,Gvarnames[i],j,Gvarnames[j]))
          }
        }
      }
      colnames(Gnet)<-c("responseNum","response","predictorNum","predictor")
    }
  }
  if(times>1){
    rlt_coef<-c()
    rlt_rs<-c()
    rlt_eq_matrix<-matrix(0,m,M)
    if(is.null(X)){
      colnames(rlt_eq_matrix)<-colnames(Y)
      rownames(rlt_eq_matrix)<-colnames(Y)
      if(!is.null(Y.fixed)){
        rlts_stability<-rlt_eq_matrix
        rlts_stability[,1:m]<-Y.fixed
      }else{
        for(s in 1:times){
          subsample<-sample(c(1:n),floor(n/2),replace=FALSE)
          Ysub<-Y[subsample,]
          ###eq_matrix saves the selected predictors for each equation
          eq_matrix<-matrix(0,m,M)
          for(k in 1:m){
            dataX=as.matrix(Ysub[,-k])
            #xy0=t(dataX)%*%Ysub[,k]
            #idx=sort.int(xy0,decreasing=TRUE,index.return=TRUE)$ix
            #dataX=dataX[,idx]
            varx<-apply(dataX,2,var)
            outl1=slim(X=dataX[,varx!=0],Y=Ysub[,k],lambda,rho,method = "lasso",verbose=FALSE)
            ##save data for linear model
            x_sel<-which(outl1$beta!=0)
            #x_sel<-idx[which(outl1$beta!=0)]
            x_num<-x_sel+as.numeric(x_sel>=k)
            eq_matrix[k,x_num]<-1  
          }
          rlt_eq_matrix<-rlt_eq_matrix+eq_matrix
        } 
        rlts_stability<-rlt_eq_matrix/times
      }
      
      ##construct equations based on stability
      if(max(rlts_stability)<stability) stop('The maximum stability we got is smaller than the stability threshold ')
      for(k in 1:m){
        if(max(rlts_stability[k,])>=stability){
          rlt_sp<-c()
          x_sel<-which(rlts_stability[k,]>=stability)
          W=as.matrix(Y[,which(rlts_stability[k,]>=stability)])
          coef=ginv(t(W)%*%W)%*%t(W)%*%Y[,k]
          sigma<-(1/n)*sum((Y[,k]-W%*%coef)^2)
          SIGMA<-sigma*ginv(t(W)%*%W)
          for(l in 1:length(x_sel)){
            Tc=coef[l]^2/diag(SIGMA)[l]
            praw<-pchisq(Tc,df=1,lower.tail=FALSE)
            #rlt_sp<-c(rlt_sp,min(praw*length(x_sel),1))##p-value for each score or variable selected in this equation
            rlt_sp<-c(rlt_sp,praw)
          }
          rlt_coef<-rbind(rlt_coef,cbind(coef,rlt_sp))
          rlt_rs<-c(rlt_rs,sigma)
          # eq_system<-paste("Y[,",k,"]~as.matrix(Y[,which(rlts_stability[",k,",]>=stability)])",sep="")
          #lm.fit<-lm(eval(parse(text=eq_system)))
          #rlt_coef<-rbind(rlt_coef,summary(lm.fit)$coefficients[-1,])
          #rlt_rs<-c(rlt_rs,sum(summary(lm.fit)$residuals^2))
        }
        else {k<-k+1} 
      }
    }
    else{
      for(s in 1:times){
        subsample<-sample(c(1:n),floor(n/2),replace=FALSE)
        Ysub<-Y[subsample,]
        Xsub<-X[subsample,]
        ##first stage
        Yhat<-matrix(0,length(subsample),m)
        xsubinv<-myinv(t(Xsub)%*%Xsub)
        Yhat<-Xsub%*%xsubinv%*%t(Xsub)%*%Ysub
        ##second stage
        ###eq_matrix saves the selected predictors for each equation
        eq_matrix<-matrix(0,m,M)
        ##beta_matrix saves the coefficients from lasso on each equation
        beta_matrix<-matrix(0,m,M-1)
        for(k in 1:m){
          dataX<-cbind(Yhat[,-k],Xsub)
          dataX<-as.matrix(dataX)
          varx<-apply(dataX,2,var)
          outl1=slim(X=dataX[,varx!=0],Y=Ysub[,k],lambda,rho,method = "lasso",verbose=FALSE)
          ##save data for linear model or linear equation
          beta_matrix[k,varx!=0]<-outl1$beta
          x_sel<-which(beta_matrix[k,]!=0)
          x_num<-x_sel+as.numeric(x_sel>=k)
          eq_matrix[k,x_num]<-1
        }
        rlt_eq_matrix<-rlt_eq_matrix+eq_matrix
      }
      ##construct equations based on stability
      rlt_p<-c()##save p-values for each gene calculated from t2 test
      rlt_sp<-c()##save p-values for each score calculated from t2 test
      YX<-cbind(Y,X)
      colnames(rlt_eq_matrix)<-colnames(YX)
      rownames(rlt_eq_matrix)<-colnames(Y)
      rlts_stability<-rlt_eq_matrix/times
      if(!is.null(Y.fixed)){
        rlts_stability[,1:m]<-Y.fixed
      }
      if(!is.null(XY.fixed)){
        rlts_stability[,(m+1):M]<-XY.fixed
      }
      if(max(rlts_stability)<stability) stop('The maximum stability we got is smaller than the stability threshold ')
      Xinv<-myinv(t(X)%*%X)
      for(k in 1:m){
        if(max(rlts_stability[k,])>=stability){
          x_sel<-which(rlts_stability[k,]>=stability)
          if(length(x_sel)>0){
            W<-YX[,x_sel]
            What<-as.matrix(X%*%Xinv%*%t(X)%*%W)
            WX<-as.matrix(t(W)%*%X%*%Xinv%*%t(X))
            coef<-ginv(WX%*%W)%*%WX%*%Y[,k]
            sigma<-(1/n)*sum((Y[,k]-W%*%coef)^2)
            SIGMA<-sigma*ginv(t(What)%*%What)
            rlt_coef<-c(rlt_coef,coef)
            rlt_rs<-c(rlt_rs,sigma)
            ##test using t2 statistics for each variable selected in this equation
            pTc<-c()
            for(l in 1:length(x_sel)){
              Tc=coef[l]^2/diag(SIGMA)[l]
              praw<-pchisq(Tc,df=1,lower.tail=FALSE)
              #rlt_sp<-c(rlt_sp,min(praw*length(x_sel),1))##p-value for each score or variable selected in this equation
              rlt_sp<-c(rlt_sp,praw)
            }
            if(sum(xsn)==1){
              rlt_p<-rlt_sp
            }else{##for cases that one gene have more than one scores
              if(is.null(names(xsn))){
                names(xsn)<-paste("x",1:length(xsn),sep="")
              }
              xsn1<-c(rep(1,m),xsn)
              xsn_num<-cumsum(xsn1)
              g_sel<-intersect(which(0<x_sel),which(x_sel<=xsn_num[1]))
              if(length(g_sel)==1){
                Tg=coef[g_sel]^2/diag(SIGMA)[g_sel]
                praw<-pchisq(Tg,df=1,lower.tail=FALSE)
                #rlt_p<-c(rlt_p,min(praw*length(x_sel),1))
                rlt_p<-c(rlt_p,praw)
              }
              if(length(g_sel)>1){
                Tg=t(coef[g_sel])%*%ginv(SIGMA[g_sel,g_sel])%*%coef[g_sel]
                praw<-pchisq(Tg,df=length(g_sel),lower.tail=FALSE)
                #rlt_p<-c(rlt_p,min(praw*length(x_sel),1))
                rlt_p<-c(rlt_p,praw)
              }
              for(j in 2:length(xsn1)){
                g_sel<-intersect(which(xsn_num[j-1]<x_sel),which(x_sel<=xsn_num[j]))
                if(length(g_sel)==1){
                  Tg=coef[g_sel]^2/diag(SIGMA)[g_sel]
                  praw<-pchisq(Tg,df=1,lower.tail=FALSE)
                  #rlt_p<-c(rlt_p,min(praw*length(x_sel),1))
                  rlt_p<-c(rlt_p,praw)
                }
                if(length(g_sel)>1){
                  Tg=t(coef[g_sel])%*%ginv(SIGMA[g_sel,g_sel])%*%coef[g_sel]
                  praw<-pchisq(Tg,df=length(g_sel),lower.tail=FALSE)
                  #rlt_p<-c(rlt_p,min(praw*length(x_sel),1))
                  rlt_p<-c(rlt_p,praw)
                }
              }
            }
          }
        }
        else{
          k=k+1
        }  
      }
    }
    
    ####save gene_score network edge
    varnames<-colnames(rlts_stability)
    Snet<-c()
    for(i in 1:nrow(rlts_stability)){
      for(j in 1:ncol(rlts_stability)){
        if(rlts_stability[i,j]>=stability){
          Snet<-rbind(Snet,c(i,varnames[i],j,varnames[j],rlts_stability[i,j]))
        }
      }
    }
    colnames(Snet)<-c("responseNum","response","predictorNum","predictor","stability")
    ###for the case that xsn is different for each gene, we should transform the the snp or score stability matrix to 
    ###gene(edge) stability matrix
    ###and based on gene stability matrix to construct gene network
    Gnet<-c()
    if(sum(xsn)==1){
      Gnet<-Snet
    }
    else{
      nx<-length(xsn)
      rlt_stability<-matrix(0,m,m+nx)
      rlt_stability[,1:m]<-rlts_stability[,1:m]
      score_num<-m
      for(i in 1:nx){
        stability_g<-as.matrix(rlts_stability[,(score_num+1):(score_num+xsn[i])])
        rlt_stability[,m+i]<-apply(stability_g,1,max)
        score_num<-score_num+xsn[i]
      }
      rownames(rlt_stability)<-colnames(Y)
      colnames(rlt_stability)<-c(colnames(Y),names(xsn))
      Gvarnames<-colnames(rlt_stability)
      for(i in 1:nrow(rlt_stability)){
        for(j in 1:ncol(rlt_stability)){
          if(rlt_stability[i,j]>=stability){
            Gnet<-rbind(Gnet,c(i,Gvarnames[i],j,Gvarnames[j],rlt_stability[i,j]))
          }
        }
      }
      colnames(Gnet)<-c("responseNum","response","predictorNum","predictor","stability")
    }
  }
  rlt<-list()
  if(times==1){
    rlt$eq_matrix<-eq_matrix
    if(sum(xsn)>1){
      rlt$G_matrix<-G_matrix
    }
  }
  else{
    rlt$eq_stability<-rlts_stability
    if(sum(xsn)>1){
      rlt$G_stability<-rlt_stability
    }
  }
  if(is.null(X)){
    rlt$eq_scorenet<-cbind(Snet,rlt_coef)
    if(times>1){
      colnames(rlt$eq_scorenet)<-c("responseNum","response","predictorNum","predictor","stability","rlt_coef","rlt_sp")
    }else{
      colnames(rlt$eq_scorenet)<-c("responseNum","response","predictorNum","predictor","rlt_coef","rlt_sp")
    }
  }
  else{
    rlt$eq_genenet<-cbind(Gnet,rlt_p)
    rlt$eq_scorenet<-cbind(Snet,rlt_coef,rlt_sp) 
  }
  
  
  rlt$eq_residuals_square<-rlt_rs
  return(rlt)
}
myinv<-function(A){
  A_svd<-fast.svd(A)
  A_inv<-A_svd$v%*%diag(1/A_svd$d)%*%t(A_svd$u)
  return(A_inv)
}

##the input of function "TableMergeByCol" are results from function "sem_effect"
##tab1 is rlt_of_sem_effect$direct
##tab2 is rlt_of_sem_effect$indirect
TableMergeByCol<-function(tab1,tab2,key1=1,key2=2){
  u1<-apply(tab1,1,function(x) paste(x[key1],x[key2],sep="\t"))
  u2<-apply(tab2,1,function(x) paste(x[key1],x[key2],sep="\t"))
  tab1<-cbind(tab1,u1)
  tab2<-cbind(tab2,u2)
  tab<-merge(tab1,tab2,by.x="u1",by.y="u2",all=T)
  tab<-tab[,c(1,4,7)]
  tab3<-as.matrix(tab[,c(2,3)])
  tab3[is.na(tab3[,1]),1]<-0
  tab3[is.na(tab3[,2]),2]<-0
  rowname<-matrix(unlist(strsplit(as.character(tab[,1]),split="\t")),ncol=2,byrow=T)
  predictor<-rowname[,1]
  response<-rowname[,2]
  rlt<-data.frame(response,predictor,tab3)
  colnames(rlt)<-c("response","predictor","direct_effect","indirect_effect")
  rlt
}

### path.effect is the result from function "sem_path_effect"
### if path.effect is not provided, the network edge information "netedge" from function "sparse_2sem" should be provided
###if the X data set is fpca score, we provide a option to give the genetic effect function for each gene
### if is.score==TRUE, the object of fpca.ksi should be provided for each fpca score
### the fpca.ksi object is the result from function "fpca.score" or function "fpca.score.geno", 
### fpca.ksi should according to the predictor fpca score
sem_effect<-function(path.effect=NULL,netedge=NULL,start=NULL,end=NULL,is.score=FALSE,fpca.score=NULL,effect.function=FALSE){
  if(is.null(path.effect)){
    if(is.null(netedge)){stop("the matrix of network edge should be provide.")}
    else{
      path.effect<-sem_path_effect(netedge,start,end)
    }
  }
  direct<-path.effect$direct
  indirect<-c()
  ind0<-c()
  for(i in 1:length(path.effect$indirect)){
    ind<-path.effect$indirect[[i]]
    if(length(ind)>0){
      for(j in 1:length(ind)){
        if(length(ind[[j]])[1]>0){
          if(is.vector(ind[[j]])){
            ind0<-rbind(ind0,ind[[j]][c(1,j+2,j+3)])
          }
          else{
            ind0<-rbind(ind0,ind[[j]][,c(1,j+2,j+3)])
          }
        } 
      }
    }
  }
  if(!is.null(ind0)){
    if(dim(ind0)[1]==1){
      indirect<-ind0
    }else{
      indirect<-ind0[!duplicated(ind0[,1:2]),]
      ind_dup<-ind0[duplicated(ind0[,1:2]),]
      if(is.vector(indirect)){
        indirect[3]<-sum(as.numeric(ind0[,3]))
      }else{
        if(is.vector(ind_dup)){
          index<-which(apply(indirect[,1:2],1,function(x) all(x==ind_dup[1:2])))
          indirect[index,3]<-as.numeric(indirect[index,3])+as.numeric(ind_dup[3])
        }else{
          for(i in 1:dim(indirect)[1]){
            index<-which(apply(ind_dup[,1:2],1,function(x) all(x==indirect[i,1:2])))
            indirect[i,3]<-as.numeric(indirect[i,3])+sum(as.numeric(ind_dup[index,3]))
          }
        }
      }
    }
  }
  colnames(direct)<-c("predictor","response","coef_direct")
  colnames(indirect)<-c("predictor","response","coef_indirect")
  if(is.score==TRUE){
    
  }
  rlt<-list()
  rlt$direct<-direct
  rlt$indirect<-indirect
  return(rlt)
  
}

###search direct and indirect effect paths from the structure equation network
###if start point isn't specified, the program will search all the paths in the network
###here the input argument "netedge" is the network matrix "eq_scorenet" from function "sparse_2sem"
sem_path_effect<-function(netedge=NULL,start=NULL,end=NULL){
  if(dim(netedge)[2]==7){
    netedge<-netedge[,which(colnames(netedge)!="stability")]
  }
  netedge<-as.matrix(netedge)
  direct<-c()
  indirect<-list()
  predictor<-unique(netedge[,"predictor"])
  if(is.null(start)){
    for(i in 1:length(predictor)){
      start<-predictor[i]
      print(paste("searching paths start from ",start))
      indirect[[i]]<-list()
      l<-1
      response1<-netedge[netedge[,"predictor"]==start,]
      if(length(response1)==dim(netedge)[2]){response1<-t(as.matrix(response1))}
      direct<-rbind(direct,response1[,c(4,2,5)])
      mideffect<-as.numeric(response1[,5])
      midpath<-cbind(start,response1[,"response"])
      while(!is.null(response1)){
        rlt_mid<-indpath(response1,mideffect,midpath,netedge,predictor)
        response1<-rlt_mid$response
        mideffect<-rlt_mid$mideffect
        midpath<-rlt_mid$midpath
        if(!is.null(response1)){
          response1<-as.matrix(response1)
          indirect[[i]][[l]]<-cbind(rlt_mid$midpath,rlt_mid$mideffect)
        }
        l<-l+1
      }
    }
    names(indirect)<-predictor
  }else{
    if(is.null(end)){
      print(paste("searching paths start from ",start))
      #indirect$predictor[i]<-list()
      l<-1
      response1<-netedge[netedge[,"predictor"]==start,]
      if(length(response1)==dim(netedge)[2]){
        response1<-t(as.matrix(response1))
      }
      direct<-rbind(direct,response1[,c(4,2,5)])
      mideffect<-as.numeric(response1[,5])
      midpath<-cbind(start,response1[,"response"])
      indirect[[1]]<-list()
      while(!is.null(response1)){
        rlt_mid<-indpath(response1,mideffect,midpath,netedge,predictor)
        response1<-rlt_mid$response
        mideffect<-rlt_mid$mideffect
        midpath<-rlt_mid$midpath
        if(!is.null(response1)){
          response1<-as.matrix(response1)
          indirect[[1]][[l]]<-cbind(rlt_mid$midpath,rlt_mid$mideffect)
        }
        l<-l+1
      }
      names(indirect)<-start
    }else{
      print(paste("searching paths start from ", start, "end to",end, "."))
      #indirect$predictor[i]<-list()
      l<-1
      response1<-netedge[netedge[,"predictor"]==start,]
      if(length(response1)==dim(netedge)[2]){
        response1<-t(as.matrix(response1))
      }
      direct<-rbind(direct,response1[response1[,"response"]==end,c(4,2,5)])
      mideffect<-as.numeric(response1[,5])
      midpath<-cbind(start,response1[,"response"])
      indirect[[1]]<-list()
      while(!is.null(response1)){
        rlt_mid<-indpath(response1,mideffect,midpath,netedge,predictor)
        response1<-rlt_mid$response
        mideffect<-rlt_mid$mideffect
        midpath<-rlt_mid$midpath
        if(!is.null(response1)){
          response1<-as.matrix(response1)
          num_end<-which(rlt_mid$midpath[,dim(rlt_mid$midpath)[2]]==end)
          if(!is.null(num_end)){
            indirect[[1]][[l]]<-cbind(rlt_mid$midpath,rlt_mid$mideffect)[num_end,]
          }
        }
        l<-l+1
      }
      names(indirect)<-start
    } 
  }
  rlt<-list()
  rlt$direct<-direct
  rlt$indirect<-indirect
  return(rlt)
}
indpath<-function(response1,mideffect,midpath,netedge,predictor){
  mideffect0<-c()
  response0<-c()
  midpath0<-c()
  for(j in 1:length(response1[,2])){
    if(response1[j,2] %in% predictor){
      response_2<-netedge[netedge[,4]==response1[j,2],]
      if(length(response_2)==dim(netedge)[2]){
        if(!response_2[2] %in% midpath[j,]){
          response2<-response_2
        }
        else{response2<-c()}
      }
      else{
        response2<-response_2[!response_2[,2]%in%midpath[j,],]
      }
      if(length(response2)>0){
        if(length(response2)==dim(netedge)[2]){
          mideffect0<-c(mideffect0,as.numeric(response2[5])*mideffect[j])
          midpath0<-rbind(midpath0,c(midpath[j,],response2[2]))
          
        }
        else{
          mideffect0<-c(mideffect0,as.numeric(response2[,5])*mideffect[j])
          midpath_0<-matrix(rep(midpath[j,],length(response2[,2])),nrow=length(response2[,2]),byrow=TRUE)
          midpath0<-rbind(midpath0,cbind(midpath_0,response2[,2]))
        }
        response0<-rbind(response0,response2)
      }
    } 
  }
  return(list("midpath"=midpath0,"mideffect"=mideffect0,"response"=response0))
}


##gene.map, file that contains gene region information with four columns: chrom, StartPos	EndPos, gene_name;
##rawdata, raw data file from plink, contain coded genotype and with first 6 columns: FID IID PAT MAT SEX PHENOTYPE
##         and then followed by snps;
##mapdata is the .map file from plink, with 4 columns: chromosome (1-22, X, Y or 0 if unplaced), rs# or snp identifier,
##         Genetic distance (morgans),Base-pair position (bp units);
##gene.extend, number of base pairs that can be extended from original gene region to be considered as one gene;
##percentage, percentage of variance should be explained by top seleted principle components;

fpca.score.geno<-function(gene.map,rawdata,mapdata,gene.extend,percentage=0.8,nbasis=37){
  xscore<-c()
  #xksi<-c()
  xprop<-c()
  xsn<-c()
  chr<-unique(gene.map[,1])
  for(i in 1:length(chr)){
    ind<-which(gene.map[,1]==chr[i])
    chrgene<-gene.map[ind,]
    snpind<-which(mapdata[,1]==chr[i])
    chrmap<-mapdata[snpind,]
    chrraw<-rawdata[,c(1:6,snpind+6)]
    snpnum<-nrow(chrmap)
    genenum<-nrow(chrgene)
    k<-1
    for(j in 1:genenum){
      x<-c()
      pos<-c()
      while((k<snpnum)&&(chrmap[k,4]<=(chrgene[j,2]-gene.extend))){k<-k+1}
      while((k<=snpnum)&&(chrmap[k,4]>(chrgene[j,2]-gene.extend))&&(chrmap[k,4]<(chrgene[j,3]+gene.extend))){
        x<-cbind(x,chrraw[,k+6])
        pos<-c(pos,chrmap[k,4])
        k<-k+1
      }
      if(!is.null(x)){
        x.matrix<-as.matrix(x)
        miss <- which(is.na(x.matrix), arr.ind = TRUE)
        item.means <- colMeans(x.matrix, na.rm = TRUE)
        x.matrix[miss] <- item.means[miss[, 2]]
        x<-x.matrix
        x_var<-apply(x,2,var)
        x<-as.matrix(x[,x_var!=0])
        pos=pos[x_var!=0]
        if(!is.null(x)){
          if(dim(x)[2]>=2){
            x_rlt<-fpca.score(x,pos=pos,gename=chrgene[j,4],percentage,nbasis)
            xscore<-cbind(xscore,x_rlt$scores)
            #xksi<-cbind(xksi,x_rlt$ksi)
            xprop<-c(xprop,x_rlt$prop)
            xsn<-c(xsn,x_rlt$xsn)
          }
        }
      }
    }
  }
  if(sum(xsn)==length(xsn)){xsn=1}
  rlt<-list()
  rlt$xscore<-xscore
  #rlt$xksi<-xksi
  rlt$xprop<-xprop
  rlt$xsn<-xsn
  return(rlt)
}



fpca.score <- function(x,pos=NULL,gename,percentage=0.8,nbasis=37){
  ninds <- dim(x)[1]
  nsnps <- dim(x)[2]
  if ( is.null(pos) ){
    pos <- (0:( nsnps-1) )/(nsnps-1)
  }else {
    idx<-order(pos)
    x<-x[,idx]
    pos<-pos[idx]
    pos<- (pos-pos[1])/(pos[nsnps]-pos[1])
  }
  
  expanded<-fourier.expansion(x,nbasis,pos)
  
  coef<-t(expanded$coef-rowMeans(expanded$coef))/sqrt(ninds) 
  pca.rlt<-prcomp(coef)
  pca.rlt$scores<-coef%*%pca.rlt$rotation
  ksi<-expanded$phi%*%pca.rlt$rotation
  v0<-diag(var(pca.rlt$scores))
  prop<-cumsum(v0)/sum(v0)
  x_sn<-as.vector(which(prop>percentage)[1])###number of scores that can explain variance larger than setted percentage
  names(x_sn)<-gename
  x_score<-pca.rlt$scores[,1:x_sn]
  x_score<-as.matrix(x_score)
  x_prop<-prop[1:x_sn]
  x_ksi<-ksi[,1:x_sn]
  x_ksi<-as.matrix(x_ksi)
  if(x_sn==1){
    colnames(x_score)<-gename
    names(x_prop)<-gename
    colnames(x_ksi)<-gename
  }
  else{
    colnames(x_score)<-paste(gename,1:x_sn,sep="_")
    names(x_prop)<-paste(gename,1:x_sn,sep="_")
    colnames(x_ksi)<-paste(gename,1:x_sn,sep="_")
  }
  rlt<-list()
  rlt$scores<-x_score
  rlt$prop<-x_prop
  rlt$ksi<-x_ksi
  rlt$xsn<-x_sn
  return(rlt)
}

fourier.expansion<- function(x,n_of_basis,pos){
  frange <- c(pos[1], pos[length(pos)])
  rlt=list();
  rlt$fbasis<-create.fourier.basis(frange,nbasis=n_of_basis)
  rlt$phi = eval.basis(pos,rlt$fbasis);
  rlt$coef<-ginv(t(rlt$phi)%*%rlt$phi)%*%t(rlt$phi)%*%t(x)
  return(rlt)
}

fourier.expansion.smoothed<- function(x,n_of_basis,pos,lambda){
  frange <- c(pos[1], pos[length(pos)])
  rlt=list();
  rlt$fbasis<-create.fourier.basis(frange,nbasis=n_of_basis)
  
  rlt$phi = eval.basis(pos,rlt$fbasis) + eval.basis(pos,rlt$fbasis,2)*lambda;
  rlt$coef<-ginv(t(rlt$phi)%*%rlt$phi)%*%t(rlt$phi)%*%t(x)
  
  return(rlt)
}

#'
#' @param adj.matrix A matrix for the structure of the network
#' @param data The normalized raw dataset
#' @param max_parent The max parent number
#' @param calculate_score The score funtion. Accept three parameters: 1.node The current node ID, 2.ParentSet The parent-set for the node and 3. The data

#sourceCpp("./src/simple_cycle.cpp")

load.score <- function( adj.matrix , data, max_parent, level, calculate_score){
  adj.graph = graph_from_adjacency_matrix(adj.matrix,mode="undirected")
  
  comps = components(adj.graph,mode="strong")
  cluster.id = comps$membership
  idx = which(comps$csize>1)
  cnns = vector(mode="list",length = length(idx) )
  cid0<-1
  for( cid in idx ){
    print(cid)
    vets = which(cluster.id==cid)
    if( length(vets) <=1 ) next;
    sg = induced.subgraph(adj.graph,vets)
    sg = is_chordal(sg,newgraph=T)$newgraph
    triangulated.mat = as_adjacency_matrix(sg)
    score.frame =  get_score(triangulated.mat, data[,vets],level[vets], calculate_score, max_parent)
    cnns[[cid0]] = list( score.frame=score.frame,vets = vets,tri.mat = triangulated.mat)
    cid0<-cid0+1
  }
  cnns
}

make.constraints <- function(triangulated.mat, score.frame ){
  
  working.parents = score.frame[,1:(ncol(score.frame)-2)]
  working.node = score.frame[,ncol(score.frame)-1]
  cons.count = nrow(score.frame)
  
  tmp.mat = as.matrix( triangulated.mat  )
  for( i in 1:nrow(tmp.mat) ){
    for( j in 1:ncol(tmp.mat) ){
      tmp.mat[i,j] = tmp.mat[i,j] * sum( working.parents[working.node==j,i] )
    }
  }
  
  tmp.mat = as.matrix( tmp.mat * 1  )
  cycles = get_cycles(tmp.mat)
  
  A = rdeque()
  
  for( cy in cycles){
    print(cy)
    cyls = cy + 1
    cyrs = c(cyls[length(cyls)],cyls[1:(length(cyls)-1)])
    
    no.pa = vector(mode="list", length=length(cyls) )
    for( i in 1:length(cyls) ){
      no.pa[[i]] = which( working.node == cyls[i] & working.parents[,cyrs[i] ] )
    }
    n=as.matrix( expand.grid(no.pa) )
    if( nrow(n) == 0) {next}
    tmp.A = matrix( 0, nrow(n),cons.count )
    for( i in 1:nrow(n) ){
      tmp.A[i, n[i,] ] = 1
    }
    A = insert_front(A, tmp.A)
  }
  return( do.call(  rbind,as.list(A) )  )
}

make.constraints2 <- function(triangulated.mat, score.frame ){
  system('rm csc.A.tmp3')
  system('touch csc.A.tmp3')

  working.parents = score.frame[,1:(ncol(score.frame)-2)]
  working.node = score.frame[,ncol(score.frame)-1]
  cons.count = nrow(score.frame)
  
  tmp.mat = as.matrix( triangulated.mat  )
  for( i in 1:nrow(tmp.mat) ){
    for( j in 1:ncol(tmp.mat) ){
      tmp.mat[i,j] = tmp.mat[i,j] * sum( working.parents[working.node==j,i] )
    }
  }
  
  tmp.mat = as.matrix( tmp.mat * 1  )
  cycles = get_cycles(tmp.mat)
  
  A = rdeque()
  
  for( cy in cycles){
    print(cy)
    cyls = cy + 1
    cyrs = c(cyls[length(cyls)],cyls[1:(length(cyls)-1)])
    
    no.pa = vector(mode="list", length=length(cyls) )
    for( i in 1:length(cyls) ){
      no.pa[[i]] = which( working.node == cyls[i] & working.parents[,cyrs[i] ] )
    }
    n=as.matrix( expand.grid(no.pa) )
    if( nrow(n) == 0) {next}
    tmp.A = matrix( 0, nrow(n),cons.count )
    for( i in 1:nrow(n) ){
      tmp.A[i, n[i,] ] = 1
    }
    A = insert_front(A, tmp.A)
  }
  for (a in as.list(A)){
    write.table(a,'csc.A.tmp3',row.names=F,col.names=F,append = T)
  }
  return(as.matrix(fread('csc.A.tmp3')))
}

make.constraints3 <- function(triangulated.mat, score.frame ){
  system('rm csc.A.tmp3')
  system('touch csc.A.tmp3')

  working.parents = score.frame[,1:(ncol(score.frame)-2)]
  working.node = score.frame[,ncol(score.frame)-1]
  cons.count = nrow(score.frame)
  
  tmp.mat = as.matrix( triangulated.mat  )
  for( i in 1:nrow(tmp.mat) ){
    for( j in 1:ncol(tmp.mat) ){
      tmp.mat[i,j] = tmp.mat[i,j] * sum( working.parents[working.node==j,i] )
    }
  }
  
  tmp.mat = as.matrix( tmp.mat * 1  )
  cycles = get_cycles(tmp.mat)
  
  A = rdeque()
  
  for( cy in cycles){
    print(cy)
    cyls = cy + 1
    cyrs = c(cyls[length(cyls)],cyls[1:(length(cyls)-1)])
    
    no.pa = vector(mode="list", length=length(cyls) )
    for( i in 1:length(cyls) ){
      no.pa[[i]] = which( working.node == cyls[i] & working.parents[,cyrs[i] ] )
    }
    n=as.matrix( expand.grid(no.pa) )
    if( nrow(n) == 0) {next}
    tmp.A = matrix( 0, nrow(n),cons.count )
    for( i in 1:nrow(n) ){
      tmp.A[i, n[i,] ] = 1
    }
    write.table(tmp.A,'csc.A.tmp3',row.names=F,col.names=F,append = T)
  }
  return(as.matrix(fread('csc.A.tmp3')))
}

get_score <- function (triangulated.mat, data,level, calculate_score, max_parent=5) {
  
  n_parms = nrow(triangulated.mat)
  index.line = n_parms+1
  score.line = n_parms+2
  
  score.frame = matrix(0,0,score.line)
  null.list = rep(0,score.line)
  
  
  for ( i in 1:n_parms ){
    if(level[i] == 0 ) next;
    possible_parent = which( (triangulated.mat[i,] > 0 ) & (level >= level[i] ) )
    
    n_p = length(possible_parent)
    n.parents = min(max_parent,n_p)
    
    
    tmp = null.list
    tmp[index.line] = i
    tmp[score.line] = sum( data[,i]^2 )
    score.frame = rbind( score.frame,tmp)
    
    if( n_p == 0){
      next
    }
    if( n_p == 1){
      tmp = null.list
      tmp[index.line] = i
      tmp[possible_parent] = 1
      tmp[score.line] = calculate_score( node=i, parent= possible_parent, data=data )
      score.frame = rbind( score.frame, tmp)
      next
    }
    
    parent.idx = vector(mode="list",length = n.parents)
    parent.scores = vector(mode="list",length = n.parents)
    for( ps in 1:n.parents){
      parents = combn(n.parents,ps)
      scores = rep(0,ncol(parents))
      parent.idx[[ps]] = lapply(seq_len(ncol(parents)),function(i) parents[,i])
      tmp = matrix(null.list,nrow=ncol(parents),ncol=score.line,byrow=T)
      keep.score = rep(T,ncol(parents))
      for( p.set in 1:ncol(parents) ){
        score = calculate_score(node=i,parent=possible_parent[parents[,p.set]],data=data)
        tmp[p.set,index.line ] = i
        tmp[p.set,possible_parent[parents[,p.set]]] = 1
        tmp[p.set,score.line] = score
        scores[p.set] = score
        if(ps > 1) {
          contains.list = unlist( lapply( parent.idx[[ps-1]], function(x) {all( x %in% parents[,p.set])} ) )
          keep.score[p.set] = all( parent.scores[[ps-1]] > score )
        }
      }
      score.frame = rbind( score.frame,tmp[keep.score,])
      parent.idx[[ps]] = parents
      parent.scores[[ps]]=scores
    }
  }
  score.frame
}

#Rcpp::sourceCpp('C:/Users/zhu2/Documents/trail/panpan/grouplasso/cnif/score_function_regression.cpp')
#Rcpp::sourceCpp('C:/Users/zhu2/Documents/trail/panpan/grouplasso/cnif/initial_sem.cpp')
#Rcpp::sourceCpp('C:/Users/zhu2/Documents/trail/panpan/grouplasso/cnif/simple_cycle.cpp')
library(igraph)

######################################
#CNIF_lower_memory
######################################

CNIF2 <- function(data,level=NULL, init.adj=NULL,max_parent =5, score_function=score_function_regression,
                 max_cycle_depth=5,lambda=0.001,rho=0.5,lowmem = FALSE){

  data = data.matrix(data)

  if( is.null(init.adj) ) {
    init.adj =  (initial_sem_nofix(data,1:ncol(data),lambda=lambda,rho=rho) != 0 )
  }

  if( is.null(level) ){
    level = rep(1,ncol(data))
  }

  colnames( init.adj ) = colnames(data)
  rownames( init.adj ) = colnames(data)

  scores.list = load.score( init.adj, data, max_parent, level,score_function )

  rlt = init.adj
  rlt[,] = 0

  for( tmp in scores.list){
    vets = tmp$vets
    score.frame = tmp$score.frame
    parms = length(vets)


    obj = score.frame[,parms+2]
    types = rep("I",parms)


    rhs   = hash()
    sense   = hash()
    A     = hash()

    working.node = score.frame[,parms+1]
    working.parent = score.frame[,1:parms]
    for( i in 1:parms ){
      cons.name = paste("one parent set",i)
      rhs[[cons.name]] = 1
      sense[[cons.name]] = "=="
      A[[cons.name]]   = (working.node == i)
    }


    cons.names = keys(A)

    A = t(values(A,cons.names))
    sense = values(sense,cons.names)
    rhs = values(rhs,cons.names)

    if(lowmem){
      write.table(A+0,'A.tmp3',row.names=F,col.names=F)
      csc.A = make.constraints3(tmp$tri.mat,score.frame )
      system('cat csc.A.tmp3 >> A.tmp3')
      A = as.matrix(fread('A.tmp3'))  
    } else {
      csc.A = make.constraints(tmp$tri.mat,score.frame )
      #A <- rbind(A,csc.A)
      mat = matrix(0,nrow=nrow(A)+nrow(csc.A),ncol=ncol(A))
      mat[1:nrow(A),] <- A
      mat[-1:-nrow(A),] <- csc.A
      A = mat
    }
    system('rm *.tmp3')
    
    sense = c(sense, rep("<=",nrow(csc.A)))
    rhs = c(rhs, rowSums(csc.A)-1)

    result = Rsymphony_solve_LP(  obj=obj,mat=A,dir=sense,rhs=rhs,types=types,verbosity = 0  )
    solution = score.frame[result$solution == 1,]

    for( i in 1:nrow(solution)){
      rlt[   vets[ solution[i,parms+1] ],
             vets[ as.logical( solution[i,1:parms]) ]
             ] = 1
    }
  }
  rlt
}