
rm(list=ls())

#methylation to methylation in clusters
load("C:/Users/zhu2/Documents/getpathway/model20170215/methylation_net/methylation_site_network_new.rda")
rlt.ms2ms <- list(data=data.grp,adj=rlt2)

#expression to expression cross all
load('C:/Users/zhu2/Documents/getpathway/model20170215/rlt_cluster_network.rda')
rlt.g2g <- list(data=exprincluster,adj=rlt_all_split,adj_c2c=cluster_network)

#qtl, methylation or snp to expression
load('C:\\Users\\zhu2\\Documents\\getpathway\\model20170215\\tacc\\qtlresult.rda')
rlt.qtlms2g <- list(rlt=rlt)
