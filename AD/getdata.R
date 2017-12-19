library(rvest)

kegg <- read_html('http://www.kegg.jp/dbget-bin/www_bget?hsa05010')
kegg <- kegg %>% html_nodes('tr:nth-child(9) .td30 div') %>% html_text()
genelist1 <- substr(kegg[1:171*2],1,regexpr(';',kegg[1:171*2])-1)

kegg <- readLines('http://www.kegg.jp/kegg-bin/show_pathway?map=hsa05010&show_description=show')
kegg <- kegg[(grep('mapdata',kegg)[2]+1):(grep('</map>',kegg)-1)]
kegg <- kegg[!grepl('onmouseover',kegg)]
kegg <- substr(kegg,regexpr('title',kegg)+7,nchar(kegg))
kegg <- substr(kegg,1,regexpr('\"',kegg)-1)
kegg <- strsplit(gsub(',','',kegg),' ')
genelist2 <- lapply(kegg,function(x){
  x <- x[1:(length(x)/2)*2]
  substr(x,2,nchar(x)-1)
})
genelist2 <- unique(genelist2)

expinpath <- do.call(cbind,lapply(expinpath,scale))
exp1 <- expinpath[,match(genelist1,colnames(expinpath))]
exp2 <- lapply(genelist2,function(x){
  x <- expinpath[,match(x,colnames(expinpath)),drop=F]
  x[is.na(x)] <- 0
  x[,apply(x,2,var)>0,drop=F]
})
exp2 <- exp2[sapply(exp2,ncol)>0]
exp2.gene <- sapply(exp2,function(x){colnames(x)[1]})
exp2.prop <- sapply(exp2,function(x){qpca(x)$prop[1]})
exp2 <- sapply(exp2,function(x) qpca(x)$X[,1,drop=F])
colnames(exp2) <- exp2.gene
save(genelist1,genelist2,exp1,exp2,file='AD_data.rda')
