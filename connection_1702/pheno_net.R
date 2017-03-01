
rm(list=ls())
setwd('C:\\Users\\zhu2\\Documents\\signaling\\')

################################
# Macro Setup
################################

# setwd('C:\\Users\\WenluluSens\\Documents\\uthealth\\signaling\\')
load('holiday.RData')
# gene31 <- paste(read.csv('clipboard')[,2])
source('codes/sparse_2sem_final.R')
source('codes/local_cnif_macro.R')
source('codes/flm_and_cca.R')
source('codes/CNIF.R')
sourceCpp("codes/score_function_regression.cpp")
sourceCpp("codes/simple_cycle.cpp")
sourceCpp("codes/initial_sem.cpp")
source('codes/CNIF_grouplasso2.R')
# source('codes/admm_group_pheno_geno_lasso.R')

plotnet <- function(x,mode='undirected'){
  diag(x) <- 0
  plot(graph_from_adjacency_matrix(t(x),mode=mode),
       edge.arrow.size=.1,
       vertex.size=3,
       vertex.label.cex=1,
       edge.width=.1)
}
det2 <- function(xcov,lambda=1){
  xcov <- xcov+diag(ncol(xcov))*lambda
  det(xcov)
}
fc <- function(x){
  w<-as.vector(t(x))[t(x)>0]
  x <- graph_from_adjacency_matrix(x>0,mode='undirected')
  fc <- membership(fastgreedy.community(x,weight=w))
  fc[] <- match(fc,unique(fc))
  fc
}
plotclust <- function(x,membership=NULL){
  G <- graph_from_adjacency_matrix(x>0)
  if(is.null(membership)){membership=rep(1,ncol(x))}
  plot(create.communities(G, membership), 
       # as.undirected(G), 
       as.directed(G),
       layout=layout.kamada.kawai(as.undirected(G)),
       edge.arrow.size=.1,
       vertex.size=3,
       vertex.label.cex=1,
       edge.width=.1)
}

###############################
# Phenotype Network
###############################

# load("C:/Users/WenluluSens/Documents/uthealth/signaling/forlinnan.rda")
varname <- c(strsplit('G_Cognitive_Score,Percept_Ori,Percept_Speed,Semantic_Memory,Working_Memory,Episodic_Memory,CHL,HDL,HDLratio,LDL,BMI,Weight,SBP,DBP,MMSE,Age,Physical_Activity,APOE,Alcohol,Smoking',',')[[1]],
             gsub(' pathway','',colnames(mpppp)[-1:-20]))
phe <- phe_data[,1:15]
sem1_phe <- sparse_2sem(phe,lambda=0.06,times=100)
plotnet(sem1_phe[[1]]>0.8)
colnames(phe) <- varname[1:ncol(phe)]
sem1_phe_adj <- CNIF(phe,init.adj=sem1_phe[[1]]>0.8,max_parent = 4)
plotnet(sem1_phe_adj)
sem1_phe_adj2 <- sem1_phe_adj;sem1_phe_adj2[rownames(sem1_phe_adj)=='chlstrl_lastvisit',colnames(sem1_phe_adj)=='bmi_lv']<-1
plotnet(sem1_phe_adj2)
set.seed(12345)
sem2_phe <- sparse_2sem(phe_data[,1:15],lambda=0.03,times=100,Y.fixed=sem1_phe_adj2,X=as.matrix(phe_data[,16:20]))
sem2_phe_adj <- rbind(sem2_phe[[1]]>0.7,0,0,0,0,0)
rownames(sem2_phe_adj) <- varname[1:20]
plotnet(sem2_phe_adj)

#################

load('C:\\Users\\zhu2\\Documents\\getpathway\\model20170215\\pid.rda')
phe2 <- phe_data_fill[match(pid,phe_data_fill$projid),match(colnames(phe_data),colnames(phe_data_fill))]

all.equal(colnames(phe2), colnames(sem2_phe_adj))
dimnames(phe2) <- list(pid,varname[1:20])
 
phenet <- list(data=phe2,adj=sem2_phe_adj)
save(phenet,file='C:\\Users\\zhu2\\Documents\\getpathway\\model20170215\\phenet.rda')
