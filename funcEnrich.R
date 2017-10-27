# Functional enrichment of genes which expression correlate with Braak regions

setwd("C:/Users/dkeo/surfdrive/Parkinson")

source("PD/base_script.R")

library("RDAVIDWebService")
# load("resources/correlated_genes.RData")
load("braakRelatedGenes.RData")

#Functional enrichment genes correlated in Braak 1-3, 4-6, or 1-6
david<-DAVIDWebService$new(email="D.L.Keo@tudelft.nl",
                           url="https://david.abcc.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
setAnnotationCategories(david, c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL"))
bg_list <- probeInfo$entrez_id
bg <- addList(david, bg_list, idType = "ENTREZ_GENE_ID", listName = "AHBA background", listType = "Background")
bg
t <- 0.05 # EASE p-value threshold
setTimeOut(david, 200000)

# Enrichment of genes correlated across braak stages
lapply(names(relatedGenes), function(b){
  genes <- relatedGenes[[b]]
  result <- addList(david, genes, idType = "ENTREZ_GENE_ID", listName = b, listType = "Gene")
  print(result)
  setCurrentBackgroundPosition(david, 1)
  getFunctionalAnnotationChartFile(david, paste0("Functional_analyses_braak/", b, "_goterms.txt"), threshold=t, count=2L)
  getClusterReportFile(david, paste0("Functional_analyses_braak/", b, "_termclusters.txt"), type = c("Term"))
})

#Read DAVID output


#Function to read Rdavid output
read.RdavidOutput <- function(fileName){
  if (file.exists(fileName)){
    terms <- read.csv(fileName, header = TRUE, sep = "\t", colClasses = "character")
    if (nrow(terms) == 0){
      print("...Removed")
      file.remove(fileName)
      NULL
    } else {
      terms
    }
  } else {
    NULL
  }
}

#Benjamini-corrected GO terms
correctedTerms <- lapply(names(relatedGenes), function(b){
  file <- paste0("Functional_analyses_braak/", b, "_goterms.txt")
  terms <- read.csv(file, header = TRUE, sep = "\t", colClasses = "character")
  rows <- which(terms$Benjamini < 0.05)
  terms <- terms[rows, c("Category", "Term", "Bonferroni", "Benjamini")]
  terms
})

#print top 3 significant terms in merged braak stages
top5terms <- lapply(correctedTerms, function(b){
  b[c(1:5),]
})
do.call(rbind, top5terms)
top5merged <- do.call(rbind, top5terms)
format(as.numeric(top5merged$Benjamini), digits = 3, scientific = TRUE)