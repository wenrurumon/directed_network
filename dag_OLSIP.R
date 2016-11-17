
#########################################################
# DAG
#########################################################

#split original graph into several non connected graphs
#provide all potential subpaths in the graph under the limitation of max parent
#provide all information necessary for the cost function calculation

subgraph <- function(adj.matrix,data,max.parent=3){
  adj.graph <- graph_from_adjacency_matrix(adj.matrix,mode='undirected')
  comps <- components(adj.graph,mode='strong')
  sub.graph <- lapply(unique(comps$membership)[comps$csize>1],function(i){
    subdata <- data[,comps$membership==i,drop=F]
    subgraph <- induced.subgraph(adj.graph,which(comps$membership==i))
    subgraph <- as_adjacency_matrix(is_chordal(subgraph,newgraph=T)$newgraph)
    subpath <- do.call(rbind,
                       lapply(1:nrow(subgraph),function(i){
                         cbind(subpath(subgraph[i,],max.parent = max.parent),i)
                       }))
    return(list(subdata=subdata,subgraph=subgraph,subpath=subpath))
  })
  return(list(info=comps,subgraph=sub.graph))
}

#provide all potential subpaths in the graph under the limitation of max parent
subpath <- function(pathi, max.parent = 3){
  n.parent <- min(sum(pathi),max.parent)
  path0 <- rep(0,length(pathi))
  if(sum(pathi)==0){
    matrix(path0,nrow=1)
  } else {
    paths <- do.call(rbind,
                     lapply(1:n.parent,function(np){
                       sel <- matrix(which(pathi>0)[combn(sum(pathi),np)],nrow=np)
                       t(sapply(1:ncol(sel),function(seli){
                         temp <- rep(0,length(pathi))
                         temp[sel[,seli]] <- 1
                         temp
                       }))
                     })
    )
    rbind(path0,paths)
  }
}

#ols cost designed for specific subgraph object
get_score_lm <- function(x.list){
  x.data <- x.list$subdata  
  x.path <- x.list$subpath
  sapply(1:nrow(x.path),function(i){
    yi <- x.data[,x.path[i,ncol(x.data)+1],drop=F]
    xi <- cbind(1,x.data[,x.path[i,1:ncol(x.data)]>0,drop=F])
    sum((lm(yi~xi-1)$residual)^2)
  })
}

#given the network framework, potential paths and cost for each path, provide directed subnetwork
ip <- function(x.graph,x.path,x.cost){
  
  parms <- ncol(x.graph)
  x.p <- data.frame(response=x.path[,ncol(x.path)],
                    nparent=rowSums(x.path[,2:ncol(x.path)-1,drop=F]),
                    cost=x.cost)
  x.p <- as.matrix(x.p %>% group_by(response,nparent) %>% summarise(min=min(cost)))
  
  keep.path <- rep(T,nrow(x.path))
  for(i in 1:nrow(x.p)){
    if(x.p[i,2]>1){
      x.sel <- (rowSums(x.path[,2:ncol(x.path)-1,drop=F])==x.p[i,2])&(x.path[,ncol(x.path)]==x.p[i,1])
      keep.path[x.sel] <- x.cost[x.sel]<x.p[i-1,3]
    }
  }
  
  x.path <- x.path[keep.path,]
  obj <- x.cost[keep.path]
  types <- rep('I',parms)
  
  rhs   = hash()
  sense   = hash()
  A     = hash()
  
  working.node = x.path[,parms+1]
  working.parent = x.path[,1:parms]
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
  csc.A = make.constraints(x.graph,cbind(x.path,obj))
  
  A = rbind(A, csc.A)
  sense = c(sense, rep("<=",nrow(csc.A)))
  rhs = c(rhs, rowSums(csc.A)-1)
  
  result = Rsymphony_solve_LP(obj=obj,mat=A,dir=sense,rhs=rhs,types=types,verbosity = 0 )
  dag <- x.path[result$solution==1,1:parms,drop=F]
  return(dag)
}


#process

x.graphs <- subgraph(adj.matrix=sem>=.8,data=Y,max.parent=5)
x.costs <- lapply(x.graphs$subgraph,get_score_lm)
x.ip <- lapply(1:length(x.graphs$subgraph),function(i){
  ip(x.graphs$subgraph[[i]]$subgraph,x.graphs$subgraph[[i]]$subpath,x.costs[[i]])
})
adj <- array(0,dim=dim(sem))
j <- 0
for(i in unique(x.graphs$info$membership)[x.graphs$info$csize>1]){
  j <- j+1
  adj[x.graphs$info$membership==i,x.graphs$info$membership==i] <- x.ip[[j]]
}
mat.sds(adj.group(adj*Y.prop,Y.group))
