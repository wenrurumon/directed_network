

##################################################################
##################################################################
##################################################################

#Score List

scores.list = load.score( rlt_sem_Y$adj>=1, expr, 3, rep(1,ncol(expr)), score_function_regression )

##################################################################
##################################################################
##################################################################

#Load Score
adj.matrix <- rlt_sem_Y$adj==1
data <- expr
max_parent <- 5
level <- rep(1,ncol(expr))
calculate_score <- score_function_regression
#Test
adj.graph <- graph_from_adjacency_matrix(adj.matrix,mode='undirected')
comps <- components(adj.graph,mode='strong')
sub.graph <- lapply(unique(comps$membership),function(i){
  subdata <- data[,comps$membership==i,drop=F]
  subgraph <- induced.subgraph(adj.graph,which(comps$membership==i))
  subgraph <- as_adjacency_matrix(is_chordal(subgraph,newgraph=T)$newgraph)
  list(subdata=subdata,subgraph=subgraph)
})
score.list <- lapply(sub.graph,function(x){
  get_score(triangulated.mat=x$subgraph,data=x$subdata,level=rep(1,ncol(x$subgraph))
            ,score_function_regression,max_parent=3)
})


##################################################################
##################################################################
##################################################################
