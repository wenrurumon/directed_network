
ref <- read.table('clipboard',header=T)
g.ref <- graph_from_data_frame(ref)
from <- sem[1,1]
to <- sem[1,2]

valii <- function(from,to,g.ref){
  ifrom <- which(names(V(g.ref))==from)
  ito <- which(names(V(g.ref))==to)
  shortest_paths(g.ref,V(g.ref)[ifrom],V(g.ref)[ito])
}
test <- lapply(1:nrow(sem),function(i){
  from <- sem[i,1]
  to <- sem[i,2]
  valii(from,to,g.ref)$vpath[[1]]
})
mean(sapply(test,length)>0)
