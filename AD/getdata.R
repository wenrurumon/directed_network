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
