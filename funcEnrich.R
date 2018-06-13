# Functional enrichment of genes which expression correlate with Braak regions

setwd("C:/Users/dkeo/surfdrive/pd_braak")

source("PD/base_script.R")
load("resources/modules.RData")
load( "resources/braakGenes.RData")
# R version should be 3.25 instead of 3.3 for RDavid
library("RDAVIDWebService")

braak <- lapply(c(positive = "pos", negative = "neg"), function(x){
  braakGenes$entrez_id[braakGenes$dir %in% x]
})

#Functional enrichment of  genes correlated greater or smaller than 0
david<-DAVIDWebService$new(email="D.L.Keo@tudelft.nl",
                           url="https://david.abcc.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
setAnnotationCategories(david, c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL"))
bg_list <- ahba.genes()
bg <- addList(david, bg_list, idType = "ENTREZ_GENE_ID", listName = "AHBA background", listType = "Background")
bg
t <- 0.05 # EASE p-value threshold
setTimeOut(david, 200000)

# Enrichment of positively and negatively correlated progression genes
lapply(names(braak), function(r){
  genes <- braak[[r]]
  result <- addList(david, genes, idType = "ENTREZ_GENE_ID", listName = r, listType = "Gene")
  print(result)
  setCurrentBackgroundPosition(david, 1)
  getFunctionalAnnotationChartFile(david, paste0("Functional_analyses/", r, "_goterms.txt"), threshold=t, count=2L)
  getClusterReportFile(david, paste0("Functional_analyses/", r, "_termclusters.txt"), type = c("Term"))
})

# Enrichment of genes in modules enriched for progression genes
lapply(names(modules[3]), function(r){
  m <- modules[[r]]
  lapply(names(m), function(l){
    genes <- m[[l]]
    result <- addList(david, genes, idType = "ENTREZ_GENE_ID", listName = paste(r, ",", l), listType = "Gene")
    print(result)
    setCurrentBackgroundPosition(david, 1)
    getFunctionalAnnotationChartFile(david, paste0("Functional_analyses/", r, "_modules/", l, "_goterms.txt"), threshold=t, count=2L)
    getClusterReportFile(david, paste0("Functional_analyses/", r, "_modules/", l, "_termclusters.txt"), type = c("Term"))
  })
})