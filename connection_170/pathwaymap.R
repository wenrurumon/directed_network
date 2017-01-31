
rm(list=ls())

url2 <- 'http://www.genome.jp/kegg/pathway.html'
html2 <- readLines(url2)
html2 <- grep('map=hsa',html2,value=T)

url2 <- paste0('http://www.genome.jp/',
  substr(html2,regexpr("/kegg-bin",html2),regexpr("&amp",html2)-1)
)

pathwaymap <- lapply(url2,function(url){
  print(url)
  html <- readLines(url)
  html <- html[which(grepl('show_pathway',html)&grepl('rect',html)&grepl('hsa',html))]
  sapply(strsplit(html,": "),function(x){
    sapply(strsplit(x[[2]],'\"'),function(x){x[[1]]})
  })
})
names(pathwaymap) <- substr(html2,regexpr('">',html2)+2,regexpr('</a>',html2)-1)

save(pathwaymap,file="C:/Users/zhu2/Documents/getpathway/model20170126/pathwaymap.rda")
